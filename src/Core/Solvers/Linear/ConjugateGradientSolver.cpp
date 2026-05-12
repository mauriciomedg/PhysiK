#include "PhysiK/Core/Solvers/Linear/ConjugateGradientSolver.h"

#include <algorithm>
#include <cmath>

namespace PhysiK
{
    namespace
    {
        constexpr float DiagonalTolerance = 1.0e-8f;
        constexpr float DenominatorTolerance = 1.0e-12f;

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

        std::vector<float> BuildScalarJacobiInverse(const SparseBlockMatrix& matrix)
        {
            const std::size_t dimension =
                static_cast<std::size_t>(std::max(0, matrix.blockCount) * 3);
            std::vector<float> inverseDiagonal(dimension, 1.0f);

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

            return inverseDiagonal;
        }

        void ApplyPreconditioner(
            const std::vector<float>& residual,
            const std::vector<float>& inverseDiagonal,
            bool useJacobiPreconditioner,
            std::vector<float>& result)
        {
            result.resize(residual.size());
            if (!useJacobiPreconditioner)
            {
                result = residual;
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
        ConjugateGradientResult result;
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

        solution.assign(dimension, 0.0f);

        std::vector<float> matrixTimesSolution;
        matrix.Multiply(solution, matrixTimesSolution);
        if (matrixTimesSolution.size() != dimension || !IsFinite(matrixTimesSolution))
        {
            solution.clear();
            return result;
        }

        std::vector<float> residual(dimension, 0.0f);
        for (std::size_t i = 0; i < dimension; ++i)
        {
            residual[i] = rhs[i] - matrixTimesSolution[i];
        }

        result.residualNorm = Norm(residual);
        const float rhsNorm = std::max(1.0f, Norm(rhs));
        const float targetResidual = std::max(0.0f, settings.tolerance) * rhsNorm;
        if (result.residualNorm <= targetResidual)
        {
            result.converged = true;
            return result;
        }

        const std::vector<float> inverseDiagonal = BuildScalarJacobiInverse(matrix);
        std::vector<float> preconditionedResidual;
        ApplyPreconditioner(
            residual,
            inverseDiagonal,
            settings.useJacobiPreconditioner,
            preconditionedResidual);
        if (!IsFinite(preconditionedResidual))
        {
            solution.clear();
            return result;
        }

        std::vector<float> direction = preconditionedResidual;
        float residualDotPreconditioned = Dot(residual, preconditionedResidual);
        if (!IsFinite(residualDotPreconditioned) || residualDotPreconditioned <= 0.0f)
        {
            solution.clear();
            return result;
        }

        std::vector<float> matrixTimesDirection;
        for (int iteration = 0; iteration < settings.maxIterations; ++iteration)
        {
            matrix.Multiply(direction, matrixTimesDirection);
            if (matrixTimesDirection.size() != dimension || !IsFinite(matrixTimesDirection))
            {
                solution.clear();
                result.converged = false;
                return result;
            }

            const float denominator = Dot(direction, matrixTimesDirection);
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

            for (std::size_t i = 0; i < dimension; ++i)
            {
                solution[i] += alpha * direction[i];
                residual[i] -= alpha * matrixTimesDirection[i];
            }

            result.iterations = iteration + 1;
            result.residualNorm = Norm(residual);
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
                residual,
                inverseDiagonal,
                settings.useJacobiPreconditioner,
                preconditionedResidual);
            if (!IsFinite(preconditionedResidual))
            {
                solution.clear();
                result.converged = false;
                return result;
            }

            const float nextResidualDotPreconditioned = Dot(residual, preconditionedResidual);
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

            for (std::size_t i = 0; i < dimension; ++i)
            {
                direction[i] = preconditionedResidual[i] + beta * direction[i];
            }

            residualDotPreconditioned = nextResidualDotPreconditioned;
        }

        return result;
    }
}
