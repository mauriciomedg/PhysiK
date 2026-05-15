#include "PhysiK/Core/Solvers/Linear/MKLLinearSolver.h"

#if defined(PHYSIK_ENABLE_MKL)

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>

#include <mkl_lapacke.h>

namespace PhysiK
{
    namespace
    {
        using Clock = std::chrono::steady_clock;

        double MillisecondsBetween(Clock::time_point start, Clock::time_point end)
        {
            return std::chrono::duration<double, std::milli>(end - start).count();
        }

        bool IsFinite(double value)
        {
            return std::isfinite(value);
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

        bool IsSquareSystem(const SparseBlockMatrix& matrix, const std::vector<float>& rhs)
        {
            const int blockCount = std::max(0, matrix.blockCount);
            return matrix.blockCount >= 0 &&
                matrix.rowStart.size() == static_cast<std::size_t>(blockCount + 1) &&
                rhs.size() == static_cast<std::size_t>(blockCount * 3) &&
                matrix.values.size() == matrix.colIndex.size() &&
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

        std::vector<double> BuildDenseRowMajorMatrix(const SparseBlockMatrix& matrix)
        {
            const int dimension = std::max(0, matrix.blockCount) * 3;
            const std::size_t denseSize =
                static_cast<std::size_t>(dimension) * static_cast<std::size_t>(dimension);
            std::vector<double> dense(denseSize, 0.0);

            for (int rowBlock = 0; rowBlock < matrix.blockCount; ++rowBlock)
            {
                const int rowBegin = matrix.rowStart[static_cast<std::size_t>(rowBlock)];
                const int rowEnd = matrix.rowStart[static_cast<std::size_t>(rowBlock + 1)];
                for (int blockIndex = rowBegin; blockIndex < rowEnd; ++blockIndex)
                {
                    const int columnBlock = matrix.colIndex[static_cast<std::size_t>(blockIndex)];
                    if (columnBlock < 0 || columnBlock >= matrix.blockCount)
                    {
                        return {};
                    }

                    const Mat3& block = matrix.values[static_cast<std::size_t>(blockIndex)];
                    for (int rowAxis = 0; rowAxis < 3; ++rowAxis)
                    {
                        const int row = rowBlock * 3 + rowAxis;
                        for (int columnAxis = 0; columnAxis < 3; ++columnAxis)
                        {
                            const float value = GetBlockValue(block, rowAxis, columnAxis);
                            if (!IsFinite(value))
                            {
                                return {};
                            }

                            const int column = columnBlock * 3 + columnAxis;
                            dense[static_cast<std::size_t>(row) *
                                static_cast<std::size_t>(dimension) +
                                static_cast<std::size_t>(column)] += value;
                        }
                    }
                }
            }

            return dense;
        }

        std::vector<double> ToDoubleVector(const std::vector<float>& values)
        {
            std::vector<double> result;
            result.reserve(values.size());
            for (float value : values)
            {
                result.push_back(static_cast<double>(value));
            }

            return result;
        }

        void ToFloatVector(const std::vector<double>& values, std::vector<float>& result)
        {
            result.clear();
            result.reserve(values.size());
            for (double value : values)
            {
                result.push_back(static_cast<float>(value));
            }
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

        double ComputeResidualNorm(
            const SparseBlockMatrix& matrix,
            const std::vector<float>& rhs,
            const std::vector<float>& solution)
        {
            std::vector<float> product;
            matrix.Multiply(solution, product);
            if (product.size() != rhs.size())
            {
                return 0.0;
            }

            std::vector<double> residual(rhs.size(), 0.0);
            for (std::size_t i = 0; i < rhs.size(); ++i)
            {
                residual[i] = static_cast<double>(rhs[i]) - static_cast<double>(product[i]);
            }

            return Norm(residual);
        }

        LinearSolveResult SolveDenseSPD(
            std::vector<double>& dense,
            const std::vector<double>& rhs,
            std::vector<double>& solution,
            const LinearSolveSettings& settings,
            double* outSolveMilliseconds = nullptr)
        {
            LinearSolveResult result;
            solution.clear();

            const int dimension = static_cast<int>(rhs.size());
            if (dimension == 0)
            {
                result.converged = true;
                return result;
            }

            if (dense.size() != static_cast<std::size_t>(dimension * dimension))
            {
                return result;
            }

            solution = rhs;
            const Clock::time_point solveStart = Clock::now();
            const int info = LAPACKE_dposv(
                LAPACK_ROW_MAJOR,
                'U',
                dimension,
                1,
                dense.data(),
                dimension,
                solution.data(),
                1);
            const Clock::time_point solveEnd = Clock::now();
            if (outSolveMilliseconds != nullptr)
            {
                *outSolveMilliseconds = MillisecondsBetween(solveStart, solveEnd);
            }

            result.iterations = 1;
            if (info != 0 || !IsFinite(solution))
            {
                solution.clear();
                result.converged = false;
                return result;
            }

            const double rhsNorm = std::max(1.0, Norm(rhs));
            const double targetResidual = std::max(0.0f, settings.tolerance) * rhsNorm;
            result.converged = true;
            result.residualNorm = static_cast<float>(targetResidual);
            return result;
        }

        void PrintMKLSolveTimings(
            const char* pathName,
            int dimension,
            int blockCount,
            std::size_t blockNonZeroCount,
            std::size_t storedScalarEntryCount,
            double conversionMilliseconds,
            double handleCreationMilliseconds,
            double solveMilliseconds,
            double handleDestructionMilliseconds,
            double totalMilliseconds)
        {
            std::printf(
                "MKLLinearSolver diagnostics [%s]\n"
                "  dimension: %d\n"
                "  block count: %d\n"
                "  block nonzeros: %zu\n"
                "  stored scalar entries: %zu\n"
                "  SparseBlockMatrix to CSR conversion time: 0.000 ms (not used by current dense path)\n"
                "  SparseBlockMatrix to dense conversion time: %.3f ms\n"
                "  MKL handle creation time: %.3f ms\n"
                "  MKL solve time: %.3f ms\n"
                "  MKL handle destruction time: %.3f ms\n"
                "  total Solve() time: %.3f ms\n",
                pathName,
                dimension,
                blockCount,
                blockNonZeroCount,
                storedScalarEntryCount,
                conversionMilliseconds,
                handleCreationMilliseconds,
                solveMilliseconds,
                handleDestructionMilliseconds,
                totalMilliseconds);
        }
    }

    LinearSolveResult MKLLinearSolver::Solve(
        const SparseBlockMatrix& matrix,
        const std::vector<float>& rhs,
        std::vector<float>& solution,
        const LinearSolveSettings& settings)
    {
        const Clock::time_point totalStart = Clock::now();
        solution.clear();
        LinearSolveResult result;
        if (!IsSquareSystem(matrix, rhs))
        {
            return result;
        }

        const int dimension = std::max(0, matrix.blockCount) * 3;
        const std::size_t blockNonZeroCount = matrix.values.size();
        const std::size_t storedScalarEntryCount = blockNonZeroCount * 9u;

        const Clock::time_point conversionStart = Clock::now();
        std::vector<double> dense = BuildDenseRowMajorMatrix(matrix);
        const Clock::time_point conversionEnd = Clock::now();
        const std::vector<double> doubleRhs = ToDoubleVector(rhs);
        std::vector<double> doubleSolution;
        double solveMilliseconds = 0.0;
        result = SolveDenseSPD(dense, doubleRhs, doubleSolution, settings, &solveMilliseconds);
        if (!result.converged)
        {
            const Clock::time_point totalEnd = Clock::now();
            PrintMKLSolveTimings(
                "SparseBlockMatrix dense LAPACKE dposv",
                dimension,
                matrix.blockCount,
                blockNonZeroCount,
                storedScalarEntryCount,
                MillisecondsBetween(conversionStart, conversionEnd),
                0.0,
                solveMilliseconds,
                0.0,
                MillisecondsBetween(totalStart, totalEnd));
            return result;
        }

        ToFloatVector(doubleSolution, solution);
        result.residualNorm = static_cast<float>(ComputeResidualNorm(matrix, rhs, solution));
        const double rhsNorm = std::max(1.0, Norm(doubleRhs));
        const double targetResidual = std::max(0.0f, settings.tolerance) * rhsNorm;
        result.converged = result.residualNorm <= static_cast<float>(targetResidual);
        if (!result.converged)
        {
            solution.clear();
        }

        const Clock::time_point totalEnd = Clock::now();
        PrintMKLSolveTimings(
            "SparseBlockMatrix dense LAPACKE dposv",
            dimension,
            matrix.blockCount,
            blockNonZeroCount,
            storedScalarEntryCount,
            MillisecondsBetween(conversionStart, conversionEnd),
            0.0,
            solveMilliseconds,
            0.0,
            MillisecondsBetween(totalStart, totalEnd));

        return result;
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

        const Clock::time_point totalStart = Clock::now();
        const Clock::time_point conversionStart = Clock::now();
        std::vector<double> dense = BuildDenseRowMajorMatrix(matrix);
        const Clock::time_point conversionEnd = Clock::now();
        double solveMilliseconds = 0.0;
        result = SolveDenseSPD(dense, rhs, solution, settings, &solveMilliseconds);

        result.residualNorm = static_cast<float>(ComputeResidualNorm(matrix, rhs, solution));
        const double rhsNorm = std::max(1.0, Norm(rhs));
        const double targetResidual = std::max(0.0f, settings.tolerance) * rhsNorm;
        result.converged = result.residualNorm <= static_cast<float>(targetResidual);
        const Clock::time_point totalEnd = Clock::now();
        PrintMKLSolveTimings(
            "CSRMatrix dense LAPACKE dposv",
            matrix.rowCount,
            matrix.rowCount,
            matrix.values.size(),
            matrix.values.size(),
            MillisecondsBetween(conversionStart, conversionEnd),
            0.0,
            solveMilliseconds,
            0.0,
            MillisecondsBetween(totalStart, totalEnd));
        return result;
    }
}

#endif
