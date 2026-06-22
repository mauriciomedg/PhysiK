#include "PhysiK/Core/Solvers/Linear/LinearSolver.h"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>

#include "PhysiK/Core/Solvers/Linear/ConjugateGradientSolver.h"

namespace PhysiK
{
    namespace
    {
        using Clock = std::chrono::steady_clock;

        double ElapsedMilliseconds(Clock::time_point start)
        {
            return std::chrono::duration<double, std::milli>(
                Clock::now() - start).count();
        }

        bool IsFinite(float value)
        {
            return std::isfinite(value);
        }

        bool IsFinite(const Vec3& value)
        {
            return IsFinite(value.x) && IsFinite(value.y) && IsFinite(value.z);
        }

        bool IsFinite(const Mat3& matrix)
        {
            return IsFinite(matrix.columns[0]) &&
                IsFinite(matrix.columns[1]) &&
                IsFinite(matrix.columns[2]);
        }

        bool BuildInversePreconditioner(
            const SparseBlockMatrix& matrix,
            bool useJacobiPreconditioner,
            std::vector<Mat3>& inversePreconditioner)
        {
            const std::size_t blockCount =
                static_cast<std::size_t>(std::max(0, matrix.blockCount));
            if (inversePreconditioner.size() != blockCount)
            {
                inversePreconditioner.resize(blockCount);
            }

            if (!useJacobiPreconditioner)
            {
                std::fill(
                    inversePreconditioner.begin(),
                    inversePreconditioner.end(),
                    Mat3::Identity());
                return true;
            }

            constexpr float DeterminantTolerance = 1.0e-8f;
            for (int block = 0; block < matrix.blockCount; ++block)
            {
                const int blockIndex = matrix.FindBlockIndex(block, block);
                if (blockIndex < 0)
                {
                    assert(false && "linear solver matrix is missing a diagonal block");
                    return false;
                }

                const Mat3& diagonalBlock =
                    matrix.values[static_cast<std::size_t>(blockIndex)];
                if (!IsFinite(diagonalBlock))
                {
                    return false;
                }

                const float determinant = Determinant(diagonalBlock);
                if (!IsFinite(determinant) ||
                    std::abs(determinant) <= DeterminantTolerance)
                {
                    return false;
                }

                const Mat3 inverseBlock = Inverse(diagonalBlock);
                if (!IsFinite(inverseBlock))
                {
                    return false;
                }

                inversePreconditioner[static_cast<std::size_t>(block)] = inverseBlock;
            }

            return true;
        }

        CurrentLinearSolver& CurrentSolver()
        {
            static CurrentLinearSolver solver;
            return solver;
        }

#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
        LinearSolverProfileData& CurrentProfile()
        {
            static LinearSolverProfileData profile;
            return profile;
        }
#endif
    }

#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
    void ResetCurrentLinearSolverProfile()
    {
        CurrentProfile() = LinearSolverProfileData{};
    }

    LinearSolverProfileData GetCurrentLinearSolverProfile()
    {
        return CurrentProfile();
    }
#endif

    LinearSolveResult CurrentLinearSolver::Solve(
        const SparseBlockMatrix& matrix,
        const std::vector<Vec3>& rhs,
        std::vector<Vec3>& solution,
        const LinearSolveSettings& settings)
    {
        LinearSolveResult result;
        Clock::time_point timerStart = Clock::now();
        if (!BuildInversePreconditioner(
                matrix,
                settings.useJacobiPreconditioner,
                inversePreconditioner))
        {
            solution.clear();
            result.preconditionerBuildMs = ElapsedMilliseconds(timerStart);
#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
            LinearSolverProfileData& profile = CurrentProfile();
            profile.totalSolveMilliseconds = result.preconditionerBuildMs;
            profile.preconditionerSetupMilliseconds =
                result.preconditionerBuildMs;
            profile.preconditionerBuildMs = result.preconditionerBuildMs;
#endif
            return result;
        }
        result.preconditionerBuildMs = ElapsedMilliseconds(timerStart);

        const int maxIterations =
            settings.maxIterations > 0
                ? settings.maxIterations
                : static_cast<int>(rhs.size()) * 3;

        timerStart = Clock::now();
        const ConjugateGradientResult cgResult =
            SolvePreconditionedConjugateGradient(
                solution,
                matrix,
                rhs,
                maxIterations,
                settings.tolerance,
                inversePreconditioner,
                residual,
                direction,
                temp);
        result.cgTotalMs = ElapsedMilliseconds(timerStart);

        result.iterations = cgResult.iterations;
        result.residualNorm = cgResult.residualNorm;
        result.converged = cgResult.converged;
        result.cgMultiplyMs = cgResult.cgMultiplyMs;
        result.cgApplyPreconditionerMs = cgResult.cgApplyPreconditionerMs;
        result.cgDotVectorOpsMs = cgResult.cgDotVectorOpsMs;
#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
        LinearSolverProfileData& profile = CurrentProfile();
        profile.totalSolveMilliseconds =
            result.preconditionerBuildMs + result.cgTotalMs;
        profile.sparseMatrixMultiplyMilliseconds = result.cgMultiplyMs;
        profile.dotProductMilliseconds = result.cgDotVectorOpsMs;
        profile.vectorUpdateMilliseconds = 0.0;
        profile.preconditionerSetupMilliseconds =
            result.preconditionerBuildMs;
        profile.preconditionerApplyMilliseconds =
            result.cgApplyPreconditionerMs;
        profile.preconditionerBuildMs = result.preconditionerBuildMs;
        profile.cgTotalMs = result.cgTotalMs;
        profile.cgMultiplyMs = result.cgMultiplyMs;
        profile.cgApplyPreconditionerMs = result.cgApplyPreconditionerMs;
        profile.cgDotVectorOpsMs = result.cgDotVectorOpsMs;
        profile.iterations = result.iterations;
        profile.residualNorm = result.residualNorm;
        profile.converged = result.converged;
#endif
        return result;
    }

    LinearSolver& GetCurrentLinearSolver()
    {
        return CurrentSolver();
    }
}
