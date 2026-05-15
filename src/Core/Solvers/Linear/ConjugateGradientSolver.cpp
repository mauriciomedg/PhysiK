#include "PhysiK/Core/Solvers/Linear/ConjugateGradientSolver.h"

#include <algorithm>
#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
#include <chrono>
#endif
#include <cmath>

namespace PhysiK
{
    namespace
    {
        constexpr float DiagonalTolerance = 1.0e-8f;
        constexpr float DenominatorTolerance = 1.0e-12f;

#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
        using Clock = std::chrono::steady_clock;

        ConjugateGradientProfileData profile;

        double MillisecondsBetween(Clock::time_point start, Clock::time_point end)
        {
            return std::chrono::duration<double, std::milli>(end - start).count();
        }

        struct ProfileScope
        {
            Clock::time_point start;
            const ConjugateGradientResult* result = nullptr;

            ~ProfileScope()
            {
                profile.totalSolveMilliseconds = MillisecondsBetween(start, Clock::now());
                if (result != nullptr)
                {
                    profile.iterations = result->iterations;
                    profile.residualNorm = result->residualNorm;
                    profile.converged = result->converged;
                }
            }
        };
#endif

        void ResizeScratchVector(std::vector<float>& values, std::size_t size)
        {
            if (values.size() != size)
            {
                values.resize(size);
            }
        }

        bool IsFinite(float value)
        {
            return std::isfinite(value);
        }

        bool IsFinite(const std::vector<float>& values)
        {
            for (float value : values)
            {
                if (!IsFinite(value))
                {
                    return false;
                }
            }

            return true;
        }

        float Dot(const std::vector<float>& a, const std::vector<float>& b)
        {
#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
            const Clock::time_point start = Clock::now();
#endif
            float result = 0.0f;
            const std::size_t count = std::min(a.size(), b.size());
            for (std::size_t i = 0; i < count; ++i)
            {
                result += a[i] * b[i];
            }

#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
            profile.dotProductMilliseconds += MillisecondsBetween(start, Clock::now());
#endif
            return result;
        }

        float Norm(const std::vector<float>& values)
        {
            const float squared = Dot(values, values);
            if (!IsFinite(squared) || squared < 0.0f)
            {
                return 0.0f;
            }

            return std::sqrt(squared);
        }

        float GetBlockValue(const Mat3& matrix, int row, int column)
        {
            const Vec3& sourceColumn = matrix.columns[column];
            if (row == 0)
            {
                return sourceColumn.x;
            }

            if (row == 1)
            {
                return sourceColumn.y;
            }

            return sourceColumn.z;
        }

        void BuildScalarJacobiInverse(
            const SparseBlockMatrix& matrix,
            std::vector<float>& inverseDiagonal)
        {
            const std::size_t dimension =
                static_cast<std::size_t>(std::max(0, matrix.blockCount) * 3);
            ResizeScratchVector(inverseDiagonal, dimension);
            std::fill(inverseDiagonal.begin(), inverseDiagonal.end(), 1.0f);

            for (int block = 0; block < matrix.blockCount; ++block)
            {
                const int blockIndex = matrix.FindBlockIndex(block, block);
                if (blockIndex < 0)
                {
                    continue;
                }

                const Mat3& diagonalBlock = matrix.values[static_cast<std::size_t>(blockIndex)];
                const std::size_t base = static_cast<std::size_t>(block * 3);
                for (int axis = 0; axis < 3; ++axis)
                {
                    const float diagonal = GetBlockValue(diagonalBlock, axis, axis);
                    if (IsFinite(diagonal) && std::abs(diagonal) > DiagonalTolerance)
                    {
                        inverseDiagonal[base + static_cast<std::size_t>(axis)] = 1.0f / diagonal;
                    }
                }
            }
        }

        void ApplyPreconditioner(
            const std::vector<float>& residual,
            const std::vector<float>& inverseDiagonal,
            bool useJacobiPreconditioner,
            std::vector<float>& result)
        {
#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
            const Clock::time_point start = Clock::now();
#endif
            ResizeScratchVector(result, residual.size());
            if (!useJacobiPreconditioner)
            {
                for (std::size_t i = 0; i < residual.size(); ++i)
                {
                    result[i] = residual[i];
                }
#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
                profile.preconditionerApplyMilliseconds += MillisecondsBetween(start, Clock::now());
#endif
                return;
            }

            for (std::size_t i = 0; i < residual.size(); ++i)
            {
                result[i] = residual[i] * inverseDiagonal[i];
            }
#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
            profile.preconditionerApplyMilliseconds += MillisecondsBetween(start, Clock::now());
#endif
        }
    }

#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
    void ResetConjugateGradientProfile()
    {
        profile = ConjugateGradientProfileData{};
    }

    ConjugateGradientProfileData GetConjugateGradientProfile()
    {
        return profile;
    }
#endif

    ConjugateGradientResult SolveConjugateGradient(
        const SparseBlockMatrix& matrix,
        const std::vector<float>& rhs,
        std::vector<float>& solution,
        const ConjugateGradientSettings& settings)
    {
        ConjugateGradientScratch scratch;
        return SolveConjugateGradient(matrix, rhs, solution, scratch, settings);
    }

    ConjugateGradientResult SolveConjugateGradient(
        const SparseBlockMatrix& matrix,
        const std::vector<float>& rhs,
        std::vector<float>& solution,
        ConjugateGradientScratch& scratch,
        const ConjugateGradientSettings& settings)
    {
#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
        ResetConjugateGradientProfile();
        const Clock::time_point totalStart = Clock::now();
#endif
        ConjugateGradientResult result;
#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
        const ProfileScope profileScope{totalStart, &result};
#endif
        const int blockCount = std::max(0, matrix.blockCount);
        const std::size_t dimension = static_cast<std::size_t>(blockCount * 3);

        solution.clear();
        if (dimension == 0)
        {
            result.converged = true;
            return result;
        }

        if (rhs.size() != dimension ||
            matrix.rowStart.size() != static_cast<std::size_t>(blockCount + 1) ||
            !IsFinite(rhs) ||
            settings.maxIterations < 0 ||
            !IsFinite(settings.tolerance))
        {
            return result;
        }

        if (solution.size() != dimension)
        {
            solution.resize(dimension);
        }
        std::fill(solution.begin(), solution.end(), 0.0f);

#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
        Clock::time_point timedStart = Clock::now();
#endif
        matrix.Multiply(solution, scratch.matrixDirection);
#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
        profile.sparseMatrixMultiplyMilliseconds += MillisecondsBetween(timedStart, Clock::now());
#endif
        if (scratch.matrixDirection.size() != dimension || !IsFinite(scratch.matrixDirection))
        {
            solution.clear();
            return result;
        }

        ResizeScratchVector(scratch.residual, dimension);
#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
        timedStart = Clock::now();
#endif
        for (std::size_t i = 0; i < dimension; ++i)
        {
            scratch.residual[i] = rhs[i] - scratch.matrixDirection[i];
        }
#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
        profile.vectorUpdateMilliseconds += MillisecondsBetween(timedStart, Clock::now());
#endif

        result.residualNorm = Norm(scratch.residual);
        const float rhsNorm = std::max(1.0f, Norm(rhs));
        const float targetResidual = std::max(0.0f, settings.tolerance) * rhsNorm;
        if (result.residualNorm <= targetResidual)
        {
            result.converged = true;
            return result;
        }

#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
        timedStart = Clock::now();
#endif
        BuildScalarJacobiInverse(matrix, scratch.inverseDiagonal);
#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
        profile.preconditionerSetupMilliseconds += MillisecondsBetween(timedStart, Clock::now());
#endif
        ApplyPreconditioner(
            scratch.residual,
            scratch.inverseDiagonal,
            settings.useJacobiPreconditioner,
            scratch.preconditionedResidual);
        if (!IsFinite(scratch.preconditionedResidual))
        {
            solution.clear();
            return result;
        }

        ResizeScratchVector(scratch.direction, dimension);
        for (std::size_t i = 0; i < dimension; ++i)
        {
            scratch.direction[i] = scratch.preconditionedResidual[i];
        }

        float residualDotPreconditioned =
            Dot(scratch.residual, scratch.preconditionedResidual);
        if (!IsFinite(residualDotPreconditioned) || residualDotPreconditioned <= 0.0f)
        {
            solution.clear();
            return result;
        }

        for (int iteration = 0; iteration < settings.maxIterations; ++iteration)
        {
#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
            timedStart = Clock::now();
#endif
            matrix.Multiply(scratch.direction, scratch.matrixDirection);
#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
            profile.sparseMatrixMultiplyMilliseconds += MillisecondsBetween(timedStart, Clock::now());
#endif
            if (scratch.matrixDirection.size() != dimension || !IsFinite(scratch.matrixDirection))
            {
                solution.clear();
                result.converged = false;
                return result;
            }

            const float denominator = Dot(scratch.direction, scratch.matrixDirection);
            if (!IsFinite(denominator) || std::abs(denominator) <= DenominatorTolerance)
            {
                solution.clear();
                result.converged = false;
                return result;
            }

            const float alpha = residualDotPreconditioned / denominator;
            if (!IsFinite(alpha))
            {
                solution.clear();
                result.converged = false;
                return result;
            }

#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
            timedStart = Clock::now();
#endif
            for (std::size_t i = 0; i < dimension; ++i)
            {
                solution[i] += alpha * scratch.direction[i];
                scratch.residual[i] -= alpha * scratch.matrixDirection[i];
            }
#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
            profile.vectorUpdateMilliseconds += MillisecondsBetween(timedStart, Clock::now());
#endif

            result.iterations = iteration + 1;
            result.residualNorm = Norm(scratch.residual);
            if (!IsFinite(result.residualNorm))
            {
                solution.clear();
                result.converged = false;
                return result;
            }

            if (result.residualNorm <= targetResidual)
            {
                result.converged = true;
                return result;
            }

            ApplyPreconditioner(
                scratch.residual,
                scratch.inverseDiagonal,
                settings.useJacobiPreconditioner,
                scratch.preconditionedResidual);
            if (!IsFinite(scratch.preconditionedResidual))
            {
                solution.clear();
                result.converged = false;
                return result;
            }

            const float nextResidualDotPreconditioned =
                Dot(scratch.residual, scratch.preconditionedResidual);
            if (!IsFinite(nextResidualDotPreconditioned) ||
                nextResidualDotPreconditioned <= 0.0f)
            {
                solution.clear();
                result.converged = false;
                return result;
            }

            const float beta = nextResidualDotPreconditioned / residualDotPreconditioned;
            if (!IsFinite(beta))
            {
                solution.clear();
                result.converged = false;
                return result;
            }

#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
            timedStart = Clock::now();
#endif
            for (std::size_t i = 0; i < dimension; ++i)
            {
                scratch.direction[i] =
                    scratch.preconditionedResidual[i] + beta * scratch.direction[i];
            }
#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
            profile.vectorUpdateMilliseconds += MillisecondsBetween(timedStart, Clock::now());
#endif

            residualDotPreconditioned = nextResidualDotPreconditioned;
        }

#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
        profile.totalSolveMilliseconds = MillisecondsBetween(totalStart, Clock::now());
        profile.iterations = result.iterations;
        profile.residualNorm = result.residualNorm;
        profile.converged = result.converged;
#endif
        return result;
    }
}
