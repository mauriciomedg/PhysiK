#include "PhysiK/Core/Solvers/Linear/ConjugateGradientSolver.h"

#include <algorithm>
#include <cmath>

namespace PhysiK
{
    namespace
    {
        constexpr float DiagonalTolerance = 1.0e-8f;
        constexpr float DenominatorTolerance = 1.0e-12f;

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
            float result = 0.0f;
            const std::size_t count = std::min(a.size(), b.size());

            for (std::size_t i = 0; i < count; ++i)
            {
                result += a[i] * b[i];
            }

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

        void EnsureDiagonalBlockIndexCache(
            const SparseBlockMatrix& matrix,
            ConjugateGradientScratch& scratch)
        {
            if (scratch.cachedBlockCount == matrix.blockCount &&
                scratch.cachedRowStart == matrix.rowStart &&
                scratch.cachedColumnIndex == matrix.colIndex)
            {
                return;
            }

            scratch.cachedBlockCount = matrix.blockCount;
            scratch.cachedRowStart = matrix.rowStart;
            scratch.cachedColumnIndex = matrix.colIndex;
            scratch.diagonalBlockIndices.assign(
                static_cast<std::size_t>(std::max(0, matrix.blockCount)),
                -1);

            for (int block = 0; block < matrix.blockCount; ++block)
            {
                const int rowBegin =
                    matrix.rowStart[static_cast<std::size_t>(block)];

                const int rowEnd =
                    matrix.rowStart[static_cast<std::size_t>(block + 1)];

                for (int blockIndex = rowBegin; blockIndex < rowEnd; ++blockIndex)
                {
                    if (matrix.colIndex[static_cast<std::size_t>(blockIndex)] == block)
                    {
                        scratch.diagonalBlockIndices[static_cast<std::size_t>(block)] =
                            blockIndex;

                        break;
                    }
                }
            }
        }

        void BuildScalarJacobiInverse(
            const SparseBlockMatrix& matrix,
            ConjugateGradientScratch& scratch)
        {
            const std::size_t dimension =
                static_cast<std::size_t>(std::max(0, matrix.blockCount) * 3);

            ResizeScratchVector(scratch.inverseDiagonal, dimension);

            std::fill(
                scratch.inverseDiagonal.begin(),
                scratch.inverseDiagonal.end(),
                1.0f);

            EnsureDiagonalBlockIndexCache(matrix, scratch);

            for (int block = 0; block < matrix.blockCount; ++block)
            {
                const int blockIndex =
                    scratch.diagonalBlockIndices[static_cast<std::size_t>(block)];

                if (blockIndex < 0)
                {
                    continue;
                }

                const Mat3& diagonalBlock =
                    matrix.values[static_cast<std::size_t>(blockIndex)];

                const std::size_t base =
                    static_cast<std::size_t>(block * 3);

                for (int axis = 0; axis < 3; ++axis)
                {
                    const float diagonal =
                        GetBlockValue(diagonalBlock, axis, axis);

                    if (IsFinite(diagonal) &&
                        std::abs(diagonal) > DiagonalTolerance)
                    {
                        scratch.inverseDiagonal[
                            base + static_cast<std::size_t>(axis)] =
                            1.0f / diagonal;
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
            ResizeScratchVector(result, residual.size());

            if (!useJacobiPreconditioner)
            {
                for (std::size_t i = 0; i < residual.size(); ++i)
                {
                    result[i] = residual[i];
                }

                return;
            }

            for (std::size_t i = 0; i < residual.size(); ++i)
            {
                result[i] = residual[i] * inverseDiagonal[i];
            }
        }
    }

    ConjugateGradientResult SolveConjugateGradient(
        const SparseBlockMatrix& matrix,
        const std::vector<float>& rhs,
        std::vector<float>& solution,
        const ConjugateGradientSettings& settings)
    {
        ConjugateGradientScratch scratch;

        return SolveConjugateGradient(
            matrix,
            rhs,
            solution,
            scratch,
            settings);
    }

    ConjugateGradientResult SolveConjugateGradient(
        const SparseBlockMatrix& matrix,
        const std::vector<float>& rhs,
        std::vector<float>& solution,
        ConjugateGradientScratch& scratch,
        const ConjugateGradientSettings& settings)
    {
        ConjugateGradientResult result;

        const int blockCount = std::max(0, matrix.blockCount);
        const std::size_t dimension =
            static_cast<std::size_t>(blockCount * 3);

        solution.clear();

        if (dimension == 0)
        {
            result.converged = true;
            return result;
        }

        if (rhs.size() != dimension ||
            matrix.rowStart.size() != static_cast<std::size_t>(blockCount + 1) ||
            !IsFinite(rhs) ||
            settings.maxIterations <= 0 ||
            !IsFinite(settings.tolerance) ||
            settings.tolerance <= 0.0f)
        {
            return result;
        }

        solution.resize(dimension);
        std::fill(solution.begin(), solution.end(), 0.0f);

        matrix.Multiply(solution, scratch.matrixDirection);

        if (scratch.matrixDirection.size() != dimension ||
            !IsFinite(scratch.matrixDirection))
        {
            solution.clear();
            return result;
        }

        ResizeScratchVector(scratch.residual, dimension);

        for (std::size_t i = 0; i < dimension; ++i)
        {
            scratch.residual[i] =
                rhs[i] - scratch.matrixDirection[i];
        }

        result.residualNorm = Norm(scratch.residual);

        const float rhsNorm = std::max(1.0f, Norm(rhs));

        const float targetResidual = settings.tolerance * rhsNorm;

        if (result.residualNorm <= targetResidual)
        {
            result.converged = true;
            return result;
        }

        BuildScalarJacobiInverse(matrix, scratch);

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
            scratch.direction[i] =
                scratch.preconditionedResidual[i];
        }

        float residualDotPreconditioned =
            Dot(
                scratch.residual,
                scratch.preconditionedResidual);

        if (!IsFinite(residualDotPreconditioned) ||
            residualDotPreconditioned <= 0.0f)
        {
            solution.clear();
            return result;
        }

        for (int iteration = 0;
            iteration < settings.maxIterations;
            ++iteration)
        {
            matrix.Multiply(
                scratch.direction,
                scratch.matrixDirection);

            if (scratch.matrixDirection.size() != dimension ||
                !IsFinite(scratch.matrixDirection))
            {
                solution.clear();
                result.converged = false;
                return result;
            }

            const float denominator =
                Dot(
                    scratch.direction,
                    scratch.matrixDirection);

            if (!IsFinite(denominator) ||
                denominator <= DenominatorTolerance)
            {
                solution.clear();
                result.converged = false;
                return result;
            }

            const float alpha =
                residualDotPreconditioned / denominator;

            if (!IsFinite(alpha))
            {
                solution.clear();
                result.converged = false;
                return result;
            }

            for (std::size_t i = 0; i < dimension; ++i)
            {
                solution[i] +=
                    alpha * scratch.direction[i];

                scratch.residual[i] -=
                    alpha * scratch.matrixDirection[i];
            }

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
                Dot(
                    scratch.residual,
                    scratch.preconditionedResidual);

            if (!IsFinite(nextResidualDotPreconditioned) ||
                nextResidualDotPreconditioned <= 0.0f)
            {
                solution.clear();
                result.converged = false;
                return result;
            }

            const float beta =
                nextResidualDotPreconditioned /
                residualDotPreconditioned;

            if (!IsFinite(beta))
            {
                solution.clear();
                result.converged = false;
                return result;
            }

            for (std::size_t i = 0; i < dimension; ++i)
            {
                scratch.direction[i] =
                    scratch.preconditionedResidual[i] +
                    beta * scratch.direction[i];
            }

            residualDotPreconditioned =
                nextResidualDotPreconditioned;
        }

        return result;
    }
}
