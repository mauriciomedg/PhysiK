#include "PhysiK/Core/Solvers/Linear/LinearSolver.h"

#include <algorithm>
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

        float GetBlockValue(const Mat3& matrix, int row, int column)
        {
            const Vec3& sourceColumn = matrix.columns[static_cast<std::size_t>(column)];
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

        void BuildInversePreconditioner(
            const SparseBlockMatrix& matrix,
            bool useJacobiPreconditioner,
            std::vector<Mat3>& inversePreconditioner)
        {
            inversePreconditioner.assign(
                static_cast<std::size_t>(std::max(0, matrix.blockCount)),
                Mat3::Identity());

            if (!useJacobiPreconditioner)
            {
                return;
            }

            constexpr float DiagonalTolerance = 1.0e-8f;
            for (int block = 0; block < matrix.blockCount; ++block)
            {
                const int blockIndex = matrix.FindBlockIndex(block, block);
                if (blockIndex < 0)
                {
                    continue;
                }

                const Mat3& diagonalBlock =
                    matrix.values[static_cast<std::size_t>(blockIndex)];
                float inverseDiagonal[3] = {1.0f, 1.0f, 1.0f};
                for (int axis = 0; axis < 3; ++axis)
                {
                    const float diagonal = GetBlockValue(diagonalBlock, axis, axis);
                    if (IsFinite(diagonal) && std::abs(diagonal) > DiagonalTolerance)
                    {
                        inverseDiagonal[axis] = 1.0f / diagonal;
                    }
                }

                inversePreconditioner[static_cast<std::size_t>(block)] =
                    Mat3::FromColumns(
                        Vec3{inverseDiagonal[0], 0.0f, 0.0f},
                        Vec3{0.0f, inverseDiagonal[1], 0.0f},
                        Vec3{0.0f, 0.0f, inverseDiagonal[2]});
            }
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
        std::vector<Mat3> inversePreconditioner;
        std::vector<Vec3> residual;
        std::vector<Vec3> direction;
        std::vector<Vec3> temp;

        BuildInversePreconditioner(
            matrix,
            settings.useJacobiPreconditioner,
            inversePreconditioner);

        const ConjugateGradientResult cgResult =
            SolvePreconditionedConjugateGradient(
                solution,
                matrix,
                rhs,
                std::max(1, settings.maxIterations),
                settings.tolerance,
                inversePreconditioner,
                residual,
                direction,
                temp);

        LinearSolveResult result;
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
