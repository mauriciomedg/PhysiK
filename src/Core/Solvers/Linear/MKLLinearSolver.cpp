#include "PhysiK/Core/Solvers/Linear/MKLLinearSolver.h"

#if defined(PHYSIK_ENABLE_MKL)

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

#include <mkl_spblas.h>

namespace PhysiK
{
    namespace
    {
        using Clock = std::chrono::steady_clock;

        constexpr float FloatDiagonalTolerance = 1.0e-8f;
        constexpr float FloatDenominatorTolerance = 1.0e-12f;
        constexpr double DoubleDiagonalTolerance = 1.0e-14;
        constexpr double DoubleDenominatorTolerance = 1.0e-24;

        struct FloatCSR
        {
            int dimension = 0;
            std::vector<MKL_INT> rowOffsets;
            std::vector<MKL_INT> columnIndices;
            std::vector<float> values;
        };

        struct DoubleCSR
        {
            int dimension = 0;
            std::vector<MKL_INT> rowOffsets;
            std::vector<MKL_INT> columnIndices;
            std::vector<double> values;
        };

        struct MKLTiming
        {
            double csrConversionMilliseconds = 0.0;
            double handleCreationMilliseconds = 0.0;
            double cgSolveMilliseconds = 0.0;
            double handleDestructionMilliseconds = 0.0;
            double totalMilliseconds = 0.0;
        };

        double MillisecondsBetween(Clock::time_point start, Clock::time_point end)
        {
            return std::chrono::duration<double, std::milli>(end - start).count();
        }

        bool IsFinite(float value)
        {
            return std::isfinite(value);
        }

        bool IsFinite(double value)
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

        float Norm(const std::vector<float>& values)
        {
            const float squared = Dot(values, values);
            if (!IsFinite(squared) || squared < 0.0f)
            {
                return 0.0f;
            }

            return std::sqrt(squared);
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

        bool IsSquareSystem(const SparseBlockMatrix& matrix, const std::vector<float>& rhs)
        {
            const int blockCount = std::max(0, matrix.blockCount);
            return matrix.blockCount >= 0 &&
                matrix.rowStart.size() == static_cast<std::size_t>(blockCount + 1) &&
                rhs.size() == static_cast<std::size_t>(blockCount * 3) &&
                matrix.values.size() == matrix.colIndex.size() &&
                IsFinite(rhs);
        }

        bool IsSquareSystem(const CSRMatrix& matrix, const std::vector<double>& rhs)
        {
            return matrix.IsValid() &&
                matrix.rowCount == matrix.colCount &&
                rhs.size() == static_cast<std::size_t>(std::max(0, matrix.rowCount)) &&
                IsFinite(rhs) &&
                IsFinite(matrix.values);
        }

        FloatCSR BuildScalarCSR(const SparseBlockMatrix& matrix)
        {
            FloatCSR csr;
            csr.dimension = std::max(0, matrix.blockCount) * 3;
            csr.rowOffsets.reserve(static_cast<std::size_t>(csr.dimension + 1));
            csr.columnIndices.reserve(matrix.values.size() * 9u);
            csr.values.reserve(matrix.values.size() * 9u);

            for (int rowBlock = 0; rowBlock < matrix.blockCount; ++rowBlock)
            {
                const int blockRowBegin = matrix.rowStart[static_cast<std::size_t>(rowBlock)];
                const int blockRowEnd = matrix.rowStart[static_cast<std::size_t>(rowBlock + 1)];
                for (int rowAxis = 0; rowAxis < 3; ++rowAxis)
                {
                    csr.rowOffsets.push_back(static_cast<MKL_INT>(csr.values.size()));
                    for (int blockIndex = blockRowBegin; blockIndex < blockRowEnd; ++blockIndex)
                    {
                        const int columnBlock = matrix.colIndex[static_cast<std::size_t>(blockIndex)];
                        if (columnBlock < 0 || columnBlock >= matrix.blockCount)
                        {
                            csr.dimension = -1;
                            return csr;
                        }

                        const Mat3& block = matrix.values[static_cast<std::size_t>(blockIndex)];
                        for (int columnAxis = 0; columnAxis < 3; ++columnAxis)
                        {
                            const float value = GetBlockValue(block, rowAxis, columnAxis);
                            if (!IsFinite(value))
                            {
                                csr.dimension = -1;
                                return csr;
                            }

                            if (value == 0.0f)
                            {
                                continue;
                            }

                            csr.columnIndices.push_back(
                                static_cast<MKL_INT>(columnBlock * 3 + columnAxis));
                            csr.values.push_back(value);
                        }
                    }
                }
            }

            csr.rowOffsets.push_back(static_cast<MKL_INT>(csr.values.size()));
            return csr;
        }

        DoubleCSR BuildScalarCSR(const CSRMatrix& matrix)
        {
            DoubleCSR csr;
            csr.dimension = std::max(0, matrix.rowCount);
            csr.rowOffsets.reserve(static_cast<std::size_t>(csr.dimension + 1));
            csr.columnIndices.reserve(matrix.values.size());
            csr.values.reserve(matrix.values.size());

            for (int row = 0; row < matrix.rowCount; ++row)
            {
                csr.rowOffsets.push_back(static_cast<MKL_INT>(csr.values.size()));
                const int rowBegin = matrix.rowOffsets[static_cast<std::size_t>(row)];
                const int rowEnd = matrix.rowOffsets[static_cast<std::size_t>(row + 1)];
                for (int valueIndex = rowBegin; valueIndex < rowEnd; ++valueIndex)
                {
                    const std::size_t index = static_cast<std::size_t>(valueIndex);
                    const double value = matrix.values[index];
                    if (!IsFinite(value))
                    {
                        csr.dimension = -1;
                        return csr;
                    }

                    if (value == 0.0)
                    {
                        continue;
                    }

                    csr.columnIndices.push_back(static_cast<MKL_INT>(matrix.columnIndices[index]));
                    csr.values.push_back(value);
                }
            }

            csr.rowOffsets.push_back(static_cast<MKL_INT>(csr.values.size()));
            return csr;
        }

        std::vector<float> BuildJacobiInverse(const FloatCSR& matrix)
        {
            std::vector<float> inverseDiagonal(
                static_cast<std::size_t>(std::max(0, matrix.dimension)),
                1.0f);

            for (int row = 0; row < matrix.dimension; ++row)
            {
                const MKL_INT rowBegin = matrix.rowOffsets[static_cast<std::size_t>(row)];
                const MKL_INT rowEnd = matrix.rowOffsets[static_cast<std::size_t>(row + 1)];
                for (MKL_INT valueIndex = rowBegin; valueIndex < rowEnd; ++valueIndex)
                {
                    const std::size_t index = static_cast<std::size_t>(valueIndex);
                    if (matrix.columnIndices[index] != row)
                    {
                        continue;
                    }

                    const float diagonal = matrix.values[index];
                    if (IsFinite(diagonal) && std::abs(diagonal) > FloatDiagonalTolerance)
                    {
                        inverseDiagonal[static_cast<std::size_t>(row)] = 1.0f / diagonal;
                    }
                    break;
                }
            }

            return inverseDiagonal;
        }

        std::vector<double> BuildJacobiInverse(const DoubleCSR& matrix)
        {
            std::vector<double> inverseDiagonal(
                static_cast<std::size_t>(std::max(0, matrix.dimension)),
                1.0);

            for (int row = 0; row < matrix.dimension; ++row)
            {
                const MKL_INT rowBegin = matrix.rowOffsets[static_cast<std::size_t>(row)];
                const MKL_INT rowEnd = matrix.rowOffsets[static_cast<std::size_t>(row + 1)];
                for (MKL_INT valueIndex = rowBegin; valueIndex < rowEnd; ++valueIndex)
                {
                    const std::size_t index = static_cast<std::size_t>(valueIndex);
                    if (matrix.columnIndices[index] != row)
                    {
                        continue;
                    }

                    const double diagonal = matrix.values[index];
                    if (IsFinite(diagonal) && std::abs(diagonal) > DoubleDiagonalTolerance)
                    {
                        inverseDiagonal[static_cast<std::size_t>(row)] = 1.0 / diagonal;
                    }
                    break;
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

        void ApplyPreconditioner(
            const std::vector<double>& residual,
            const std::vector<double>& inverseDiagonal,
            bool useJacobiPreconditioner,
            std::vector<double>& result)
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

        bool Multiply(sparse_matrix_t handle, const std::vector<float>& input, std::vector<float>& output)
        {
            output.assign(input.size(), 0.0f);
            matrix_descr descriptor;
            descriptor.type = SPARSE_MATRIX_TYPE_GENERAL;
            descriptor.mode = SPARSE_FILL_MODE_FULL;
            descriptor.diag = SPARSE_DIAG_NON_UNIT;

            const sparse_status_t status = mkl_sparse_s_mv(
                SPARSE_OPERATION_NON_TRANSPOSE,
                1.0f,
                handle,
                descriptor,
                input.data(),
                0.0f,
                output.data());
            return status == SPARSE_STATUS_SUCCESS && IsFinite(output);
        }

        bool Multiply(sparse_matrix_t handle, const std::vector<double>& input, std::vector<double>& output)
        {
            output.assign(input.size(), 0.0);
            matrix_descr descriptor;
            descriptor.type = SPARSE_MATRIX_TYPE_GENERAL;
            descriptor.mode = SPARSE_FILL_MODE_FULL;
            descriptor.diag = SPARSE_DIAG_NON_UNIT;

            const sparse_status_t status = mkl_sparse_d_mv(
                SPARSE_OPERATION_NON_TRANSPOSE,
                1.0,
                handle,
                descriptor,
                input.data(),
                0.0,
                output.data());
            return status == SPARSE_STATUS_SUCCESS && IsFinite(output);
        }

        bool CreateSparseHandle(FloatCSR& matrix, sparse_matrix_t& handle)
        {
            if (matrix.dimension <= 0 ||
                matrix.rowOffsets.size() != static_cast<std::size_t>(matrix.dimension + 1) ||
                matrix.columnIndices.size() != matrix.values.size())
            {
                return false;
            }

            const sparse_status_t status = mkl_sparse_s_create_csr(
                &handle,
                SPARSE_INDEX_BASE_ZERO,
                static_cast<MKL_INT>(matrix.dimension),
                static_cast<MKL_INT>(matrix.dimension),
                matrix.rowOffsets.data(),
                matrix.rowOffsets.data() + 1,
                matrix.columnIndices.data(),
                matrix.values.data());
            return status == SPARSE_STATUS_SUCCESS;
        }

        bool CreateSparseHandle(DoubleCSR& matrix, sparse_matrix_t& handle)
        {
            if (matrix.dimension <= 0 ||
                matrix.rowOffsets.size() != static_cast<std::size_t>(matrix.dimension + 1) ||
                matrix.columnIndices.size() != matrix.values.size())
            {
                return false;
            }

            const sparse_status_t status = mkl_sparse_d_create_csr(
                &handle,
                SPARSE_INDEX_BASE_ZERO,
                static_cast<MKL_INT>(matrix.dimension),
                static_cast<MKL_INT>(matrix.dimension),
                matrix.rowOffsets.data(),
                matrix.rowOffsets.data() + 1,
                matrix.columnIndices.data(),
                matrix.values.data());
            return status == SPARSE_STATUS_SUCCESS;
        }

        LinearSolveResult SolveCG(
            sparse_matrix_t handle,
            const FloatCSR& matrix,
            const std::vector<float>& rhs,
            std::vector<float>& solution,
            const LinearSolveSettings& settings)
        {
            LinearSolveResult result;
            const std::size_t dimension = static_cast<std::size_t>(std::max(0, matrix.dimension));
            solution.clear();
            if (dimension == 0)
            {
                result.converged = true;
                return result;
            }

            if (rhs.size() != dimension ||
                settings.maxIterations < 0 ||
                !IsFinite(settings.tolerance) ||
                !IsFinite(rhs))
            {
                return result;
            }

            solution.assign(dimension, 0.0f);
            std::vector<float> matrixTimesSolution;
            if (!Multiply(handle, solution, matrixTimesSolution))
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

            const std::vector<float> inverseDiagonal = BuildJacobiInverse(matrix);
            std::vector<float> preconditionedResidual;
            ApplyPreconditioner(residual, inverseDiagonal, settings.useJacobiPreconditioner, preconditionedResidual);
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
                if (!Multiply(handle, direction, matrixTimesDirection))
                {
                    solution.clear();
                    result.converged = false;
                    return result;
                }

                const float denominator = Dot(direction, matrixTimesDirection);
                if (!IsFinite(denominator) || std::abs(denominator) <= FloatDenominatorTolerance)
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

                ApplyPreconditioner(residual, inverseDiagonal, settings.useJacobiPreconditioner, preconditionedResidual);
                if (!IsFinite(preconditionedResidual))
                {
                    solution.clear();
                    result.converged = false;
                    return result;
                }

                const float nextResidualDotPreconditioned = Dot(residual, preconditionedResidual);
                if (!IsFinite(nextResidualDotPreconditioned) || nextResidualDotPreconditioned <= 0.0f)
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

        LinearSolveResult SolveCG(
            sparse_matrix_t handle,
            const DoubleCSR& matrix,
            const std::vector<double>& rhs,
            std::vector<double>& solution,
            const LinearSolveSettings& settings)
        {
            LinearSolveResult result;
            const std::size_t dimension = static_cast<std::size_t>(std::max(0, matrix.dimension));
            solution.clear();
            if (dimension == 0)
            {
                result.converged = true;
                return result;
            }

            if (rhs.size() != dimension ||
                settings.maxIterations < 0 ||
                !IsFinite(settings.tolerance) ||
                !IsFinite(rhs))
            {
                return result;
            }

            solution.assign(dimension, 0.0);
            std::vector<double> matrixTimesSolution;
            if (!Multiply(handle, solution, matrixTimesSolution))
            {
                solution.clear();
                return result;
            }

            std::vector<double> residual(dimension, 0.0);
            for (std::size_t i = 0; i < dimension; ++i)
            {
                residual[i] = rhs[i] - matrixTimesSolution[i];
            }

            double residualNorm = Norm(residual);
            result.residualNorm = static_cast<float>(residualNorm);
            const double rhsNorm = std::max(1.0, Norm(rhs));
            const double targetResidual = std::max(0.0f, settings.tolerance) * rhsNorm;
            if (residualNorm <= targetResidual)
            {
                result.converged = true;
                return result;
            }

            const std::vector<double> inverseDiagonal = BuildJacobiInverse(matrix);
            std::vector<double> preconditionedResidual;
            ApplyPreconditioner(residual, inverseDiagonal, settings.useJacobiPreconditioner, preconditionedResidual);
            if (!IsFinite(preconditionedResidual))
            {
                solution.clear();
                return result;
            }

            std::vector<double> direction = preconditionedResidual;
            double residualDotPreconditioned = Dot(residual, preconditionedResidual);
            if (!IsFinite(residualDotPreconditioned) || residualDotPreconditioned <= 0.0)
            {
                solution.clear();
                return result;
            }

            std::vector<double> matrixTimesDirection;
            for (int iteration = 0; iteration < settings.maxIterations; ++iteration)
            {
                if (!Multiply(handle, direction, matrixTimesDirection))
                {
                    solution.clear();
                    result.converged = false;
                    return result;
                }

                const double denominator = Dot(direction, matrixTimesDirection);
                if (!IsFinite(denominator) || std::abs(denominator) <= DoubleDenominatorTolerance)
                {
                    solution.clear();
                    result.converged = false;
                    return result;
                }

                const double alpha = residualDotPreconditioned / denominator;
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
                residualNorm = Norm(residual);
                result.residualNorm = static_cast<float>(residualNorm);
                if (!IsFinite(residualNorm))
                {
                    solution.clear();
                    result.converged = false;
                    return result;
                }

                if (residualNorm <= targetResidual)
                {
                    result.converged = true;
                    return result;
                }

                ApplyPreconditioner(residual, inverseDiagonal, settings.useJacobiPreconditioner, preconditionedResidual);
                if (!IsFinite(preconditionedResidual))
                {
                    solution.clear();
                    result.converged = false;
                    return result;
                }

                const double nextResidualDotPreconditioned = Dot(residual, preconditionedResidual);
                if (!IsFinite(nextResidualDotPreconditioned) || nextResidualDotPreconditioned <= 0.0)
                {
                    solution.clear();
                    result.converged = false;
                    return result;
                }

                const double beta = nextResidualDotPreconditioned / residualDotPreconditioned;
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

        void PrintDiagnostics(
            const char* pathName,
            int dimension,
            std::size_t nonZeroCount,
            const LinearSolveResult& result,
            const MKLTiming& timing)
        {
            std::printf(
                "MKLLinearSolver diagnostics [%s]\n"
                "  dimension: %d\n"
                "  CSR nonzeros: %zu\n"
                "  CSR conversion time: %.3f ms\n"
                "  MKL sparse handle creation time: %.3f ms\n"
                "  MKL CG solve time: %.3f ms\n"
                "  MKL sparse handle destruction time: %.3f ms\n"
                "  iterations: %d\n"
                "  residual norm: %.8g\n"
                "  total Solve() time: %.3f ms\n",
                pathName,
                dimension,
                nonZeroCount,
                timing.csrConversionMilliseconds,
                timing.handleCreationMilliseconds,
                timing.cgSolveMilliseconds,
                timing.handleDestructionMilliseconds,
                result.iterations,
                result.residualNorm,
                timing.totalMilliseconds);
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
        MKLTiming timing;
        if (!IsSquareSystem(matrix, rhs))
        {
            return result;
        }

        const Clock::time_point conversionStart = Clock::now();
        FloatCSR csr = BuildScalarCSR(matrix);
        const Clock::time_point conversionEnd = Clock::now();
        timing.csrConversionMilliseconds = MillisecondsBetween(conversionStart, conversionEnd);
        if (csr.dimension == 0)
        {
            result.converged = true;
            timing.totalMilliseconds = MillisecondsBetween(totalStart, Clock::now());
            PrintDiagnostics("SparseBlockMatrix MKL sparse CG", csr.dimension, csr.values.size(), result, timing);
            return result;
        }

        sparse_matrix_t handle = nullptr;
        const Clock::time_point handleStart = Clock::now();
        const bool handleCreated = CreateSparseHandle(csr, handle);
        const Clock::time_point handleEnd = Clock::now();
        timing.handleCreationMilliseconds = MillisecondsBetween(handleStart, handleEnd);
        if (!handleCreated)
        {
            timing.totalMilliseconds = MillisecondsBetween(totalStart, Clock::now());
            PrintDiagnostics("SparseBlockMatrix MKL sparse CG", csr.dimension, csr.values.size(), result, timing);
            return result;
        }

        const Clock::time_point solveStart = Clock::now();
        result = SolveCG(handle, csr, rhs, solution, settings);
        const Clock::time_point solveEnd = Clock::now();
        timing.cgSolveMilliseconds = MillisecondsBetween(solveStart, solveEnd);

        const Clock::time_point destroyStart = Clock::now();
        mkl_sparse_destroy(handle);
        const Clock::time_point destroyEnd = Clock::now();
        timing.handleDestructionMilliseconds = MillisecondsBetween(destroyStart, destroyEnd);
        timing.totalMilliseconds = MillisecondsBetween(totalStart, destroyEnd);

        PrintDiagnostics("SparseBlockMatrix MKL sparse CG", csr.dimension, csr.values.size(), result, timing);
        return result;
    }

    LinearSolveResult MKLLinearSolver::SolveSPD(
        const CSRMatrix& matrix,
        const std::vector<double>& rhs,
        std::vector<double>& solution,
        const LinearSolveSettings& settings)
    {
        const Clock::time_point totalStart = Clock::now();
        solution.clear();
        LinearSolveResult result;
        MKLTiming timing;
        if (!IsSquareSystem(matrix, rhs))
        {
            return result;
        }

        const Clock::time_point conversionStart = Clock::now();
        DoubleCSR csr = BuildScalarCSR(matrix);
        const Clock::time_point conversionEnd = Clock::now();
        timing.csrConversionMilliseconds = MillisecondsBetween(conversionStart, conversionEnd);
        if (csr.dimension == 0)
        {
            result.converged = true;
            timing.totalMilliseconds = MillisecondsBetween(totalStart, Clock::now());
            PrintDiagnostics("CSRMatrix MKL sparse CG", csr.dimension, csr.values.size(), result, timing);
            return result;
        }

        sparse_matrix_t handle = nullptr;
        const Clock::time_point handleStart = Clock::now();
        const bool handleCreated = CreateSparseHandle(csr, handle);
        const Clock::time_point handleEnd = Clock::now();
        timing.handleCreationMilliseconds = MillisecondsBetween(handleStart, handleEnd);
        if (!handleCreated)
        {
            timing.totalMilliseconds = MillisecondsBetween(totalStart, Clock::now());
            PrintDiagnostics("CSRMatrix MKL sparse CG", csr.dimension, csr.values.size(), result, timing);
            return result;
        }

        const Clock::time_point solveStart = Clock::now();
        result = SolveCG(handle, csr, rhs, solution, settings);
        const Clock::time_point solveEnd = Clock::now();
        timing.cgSolveMilliseconds = MillisecondsBetween(solveStart, solveEnd);

        const Clock::time_point destroyStart = Clock::now();
        mkl_sparse_destroy(handle);
        const Clock::time_point destroyEnd = Clock::now();
        timing.handleDestructionMilliseconds = MillisecondsBetween(destroyStart, destroyEnd);
        timing.totalMilliseconds = MillisecondsBetween(totalStart, destroyEnd);

        PrintDiagnostics("CSRMatrix MKL sparse CG", csr.dimension, csr.values.size(), result, timing);
        return result;
    }
}

#endif
