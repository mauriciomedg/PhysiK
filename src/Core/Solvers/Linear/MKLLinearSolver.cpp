#include "PhysiK/Core/Solvers/Linear/MKLLinearSolver.h"

#if defined(PHYSIK_ENABLE_MKL)

#include <algorithm>
#include <cmath>

#include <mkl_lapacke.h>

namespace PhysiK
{
    namespace
    {
        bool IsFinite(double value)
        {
            return std::isfinite(value);
        }

        bool IsFinite(const std::vector<double>& values)
        {
            for (double value : values)
            {
                if (!IsFinite(value))
                {
                    return false;
                }
            }

            return true;
        }

        double Dot(const std::vector<double>& a, const std::vector<double>& b)
        {
            double result = 0.0;
            const std::size_t count = std::min(a.size(), b.size());
            for (std::size_t i = 0; i < count; ++i)
            {
                result += a[i] * b[i];
            }

            return result;
        }

        double Norm(const std::vector<double>& values)
        {
            const double squared = Dot(values, values);
            if (!IsFinite(squared) || squared < 0.0)
            {
                return 0.0;
            }

            return std::sqrt(squared);
        }

        bool IsSquareSystem(const CSRMatrix& matrix, const std::vector<double>& rhs)
        {
            return matrix.IsValid() &&
                matrix.rowCount == matrix.colCount &&
                rhs.size() == static_cast<std::size_t>(std::max(0, matrix.rowCount)) &&
                IsFinite(rhs);
        }

        std::vector<double> BuildDenseRowMajorMatrix(const CSRMatrix& matrix)
        {
            const int dimension = std::max(0, matrix.rowCount);
            const std::size_t denseSize =
                static_cast<std::size_t>(dimension) * static_cast<std::size_t>(dimension);
            std::vector<double> dense(denseSize, 0.0);

            for (int row = 0; row < dimension; ++row)
            {
                const int rowBegin = matrix.rowOffsets[static_cast<std::size_t>(row)];
                const int rowEnd = matrix.rowOffsets[static_cast<std::size_t>(row + 1)];
                for (int valueIndex = rowBegin; valueIndex < rowEnd; ++valueIndex)
                {
                    const std::size_t index = static_cast<std::size_t>(valueIndex);
                    const int column = matrix.columnIndices[index];
                    dense[static_cast<std::size_t>(row) *
                        static_cast<std::size_t>(dimension) +
                        static_cast<std::size_t>(column)] += matrix.values[index];
                }
            }

            return dense;
        }

        double ComputeResidualNorm(
            const CSRMatrix& matrix,
            const std::vector<double>& rhs,
            const std::vector<double>& solution)
        {
            std::vector<double> product;
            matrix.Multiply(solution, product);
            if (product.size() != rhs.size())
            {
                return 0.0;
            }

            std::vector<double> residual(rhs.size(), 0.0);
            for (std::size_t i = 0; i < rhs.size(); ++i)
            {
                residual[i] = rhs[i] - product[i];
            }

            return Norm(residual);
        }
    }

    LinearSolveResult MKLLinearSolver::Solve(
        const SparseBlockMatrix&,
        const std::vector<float>&,
        std::vector<float>& solution,
        const LinearSolveSettings&)
    {
        solution.clear();
        return LinearSolveResult{};
    }

    LinearSolveResult MKLLinearSolver::SolveSPD(
        const CSRMatrix& matrix,
        const std::vector<double>& rhs,
        std::vector<double>& solution,
        const LinearSolveSettings& settings)
    {
        LinearSolveResult result;
        solution.clear();

        if (!IsSquareSystem(matrix, rhs) || !IsFinite(matrix.values))
        {
            return result;
        }

        const int dimension = matrix.rowCount;
        if (dimension == 0)
        {
            result.converged = true;
            return result;
        }

        std::vector<double> dense = BuildDenseRowMajorMatrix(matrix);
        solution = rhs;

        const int info = LAPACKE_dposv(
            LAPACK_ROW_MAJOR,
            'U',
            dimension,
            1,
            dense.data(),
            dimension,
            solution.data(),
            1);

        result.iterations = 1;
        if (info != 0 || !IsFinite(solution))
        {
            solution.clear();
            result.converged = false;
            return result;
        }

        result.residualNorm = static_cast<float>(ComputeResidualNorm(matrix, rhs, solution));
        const double rhsNorm = std::max(1.0, Norm(rhs));
        const double targetResidual = std::max(0.0f, settings.tolerance) * rhsNorm;
        result.converged = result.residualNorm <= static_cast<float>(targetResidual);
        return result;
    }
}

#endif
