#include "PhysiK/Core/Solvers/Linear/DenseLinearSolver.h"

#include <algorithm>
#include <cmath>

namespace PhysiK
{
    namespace
    {
        bool IsFinite(const std::vector<float>& values)
        {
            for (float value : values)
            {
                if (!std::isfinite(value))
                {
                    return false;
                }
            }

            return true;
        }
    }

    bool DenseLinearSolver::Solve(
        const std::vector<float>& matrix,
        const std::vector<float>& rhs,
        std::vector<float>& solution,
        int dimension)
    {
        solution.clear();

        if (dimension < 0)
        {
            return false;
        }

        if (dimension == 0)
        {
            return true;
        }

        const std::size_t n = static_cast<std::size_t>(dimension);
        if (matrix.size() != n * n || rhs.size() != n || !IsFinite(matrix) || !IsFinite(rhs))
        {
            return false;
        }

        std::vector<float> a = matrix;
        std::vector<float> b = rhs;

        constexpr float PivotTolerance = 1.0e-8f;

        for (int column = 0; column < dimension; ++column)
        {
            int pivotRow = column;
            float pivotAbs = std::abs(a[static_cast<std::size_t>(column) * n + column]);

            for (int row = column + 1; row < dimension; ++row)
            {
                const float candidate =
                    std::abs(a[static_cast<std::size_t>(row) * n + column]);
                if (candidate > pivotAbs)
                {
                    pivotAbs = candidate;
                    pivotRow = row;
                }
            }

            if (pivotAbs <= PivotTolerance)
            {
                return false;
            }

            if (pivotRow != column)
            {
                for (int c = column; c < dimension; ++c)
                {
                    std::swap(
                        a[static_cast<std::size_t>(column) * n + c],
                        a[static_cast<std::size_t>(pivotRow) * n + c]);
                }
                std::swap(b[static_cast<std::size_t>(column)], b[static_cast<std::size_t>(pivotRow)]);
            }

            const float pivot = a[static_cast<std::size_t>(column) * n + column];
            for (int row = column + 1; row < dimension; ++row)
            {
                const float factor = a[static_cast<std::size_t>(row) * n + column] / pivot;
                if (factor == 0.0f)
                {
                    continue;
                }

                a[static_cast<std::size_t>(row) * n + column] = 0.0f;
                for (int c = column + 1; c < dimension; ++c)
                {
                    a[static_cast<std::size_t>(row) * n + c] -=
                        factor * a[static_cast<std::size_t>(column) * n + c];
                }
                b[static_cast<std::size_t>(row)] -= factor * b[static_cast<std::size_t>(column)];
            }
        }

        solution.assign(n, 0.0f);
        for (int row = dimension - 1; row >= 0; --row)
        {
            float sum = b[static_cast<std::size_t>(row)];
            for (int column = row + 1; column < dimension; ++column)
            {
                sum -= a[static_cast<std::size_t>(row) * n + column] *
                    solution[static_cast<std::size_t>(column)];
            }

            const float diagonal = a[static_cast<std::size_t>(row) * n + row];
            if (std::abs(diagonal) <= PivotTolerance)
            {
                solution.clear();
                return false;
            }

            solution[static_cast<std::size_t>(row)] = sum / diagonal;
        }

        if (!IsFinite(solution))
        {
            solution.clear();
            return false;
        }

        return true;
    }
}
