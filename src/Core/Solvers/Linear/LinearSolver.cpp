#include "PhysiK/Core/Solvers/Linear/LinearSolver.h"

#include <algorithm>
#include <cassert>
#include <cmath>

#include "PhysiK/Core/Solvers/Linear/ConjugateGradientSolver.h"

namespace PhysiK
{
    namespace
    {
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
    }

    LinearSolveResult CurrentLinearSolver::Solve(
        const SparseBlockMatrix& matrix,
        const std::vector<Vec3>& rhs,
        std::vector<Vec3>& solution,
        const LinearSolveSettings& settings)
    {
        LinearSolveResult result;
        const int maxIterations =
            settings.maxIterations > 0
                ? settings.maxIterations
                : static_cast<int>(rhs.size()) * 3;

        if (!BuildInversePreconditioner(
                matrix,
                settings.useJacobiPreconditioner,
                inversePreconditioner))
        {
            solution.clear();
            return result;
        }

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

        result.iterations = cgResult.iterations;
        result.residualNorm = cgResult.residualNorm;
        result.converged = cgResult.converged;
        return result;
    }

    LinearSolver& GetCurrentLinearSolver()
    {
        return CurrentSolver();
    }
}
