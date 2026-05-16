#include "PhysiK/Math/SparseBlockMatrix.h"

#include <algorithm>
#if defined(PHYSIK_ENABLE_PERF_LOGGING)
#include <chrono>
#endif

namespace PhysiK
{
    namespace
    {
#if defined(PHYSIK_ENABLE_PERF_LOGGING)
        thread_local bool timingEnabled = false;
        thread_local double multiplyMilliseconds = 0.0;

        using Clock = std::chrono::steady_clock;

        double ElapsedMilliseconds(Clock::time_point start)
        {
            return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
        }
#endif

        Mat3 Add(const Mat3& a, const Mat3& b)
        {
            return Mat3::FromColumns(
                a.columns[0] + b.columns[0],
                a.columns[1] + b.columns[1],
                a.columns[2] + b.columns[2]);
        }

        bool IsValidCoordinate(int rowBlock, int colBlock, int blockCount)
        {
            return rowBlock >= 0 && rowBlock < blockCount &&
                colBlock >= 0 && colBlock < blockCount;
        }
    }

    void SetSparseBlockMatrixTimingEnabled(bool enabled)
    {
#if defined(PHYSIK_ENABLE_PERF_LOGGING)
        timingEnabled = enabled;
#else
        (void)enabled;
#endif
    }

    void ResetSparseBlockMatrixTiming()
    {
#if defined(PHYSIK_ENABLE_PERF_LOGGING)
        multiplyMilliseconds = 0.0;
#endif
    }

    double GetSparseBlockMatrixMultiplyMilliseconds()
    {
#if defined(PHYSIK_ENABLE_PERF_LOGGING)
        return multiplyMilliseconds;
#else
        return 0.0;
#endif
    }

    void SparseBlockMatrix::Clear()
    {
        blockCount = 0;
        rowStart.clear();
        colIndex.clear();
        values.clear();
        blockLookup.clear();
    }

    void SparseBlockMatrix::ClearValues()
    {
        std::fill(values.begin(), values.end(), Mat3::Zero());
    }

    void SparseBlockMatrix::BuildPattern(
        int newBlockCount,
        const std::vector<std::pair<int, int>>& blockCoordinates)
    {
        Clear();
        if (newBlockCount <= 0)
        {
            return;
        }

        blockCount = newBlockCount;
        std::vector<std::vector<int>> rows(static_cast<std::size_t>(blockCount));
        for (const std::pair<int, int>& coordinate : blockCoordinates)
        {
            const int rowBlock = coordinate.first;
            const int colBlock = coordinate.second;
            if (!IsValidCoordinate(rowBlock, colBlock, blockCount))
            {
                continue;
            }

            rows[static_cast<std::size_t>(rowBlock)].push_back(colBlock);
        }

        rowStart.assign(static_cast<std::size_t>(blockCount + 1), 0);
        for (int row = 0; row < blockCount; ++row)
        {
            std::vector<int>& columns = rows[static_cast<std::size_t>(row)];
            std::sort(columns.begin(), columns.end());
            columns.erase(std::unique(columns.begin(), columns.end()), columns.end());
            rowStart[static_cast<std::size_t>(row + 1)] =
                rowStart[static_cast<std::size_t>(row)] + static_cast<int>(columns.size());

            for (int column : columns)
            {
                const int blockIndex = static_cast<int>(colIndex.size());
                colIndex.push_back(column);
                values.push_back(Mat3::Zero());
                blockLookup[MakeKey(row, column)] = blockIndex;
            }
        }
    }

    bool SparseBlockMatrix::AddBlock(int rowBlock, int colBlock, const Mat3& block)
    {
        const int blockIndex = FindBlockIndex(rowBlock, colBlock);
        if (blockIndex < 0)
        {
            return false;
        }

        values[static_cast<std::size_t>(blockIndex)] =
            Add(values[static_cast<std::size_t>(blockIndex)], block);
        return true;
    }

    void SparseBlockMatrix::Multiply(
        const std::vector<float>& input,
        std::vector<float>& output) const
    {
#if defined(PHYSIK_ENABLE_PERF_LOGGING)
        const bool recordTiming = timingEnabled;
        const Clock::time_point start = recordTiming ? Clock::now() : Clock::time_point{};
#endif
        const std::size_t dimension = static_cast<std::size_t>(std::max(0, blockCount) * 3);
        if (input.size() < dimension || rowStart.size() != static_cast<std::size_t>(blockCount + 1))
        {
            output.assign(dimension, 0.0f);
#if defined(PHYSIK_ENABLE_PERF_LOGGING)
            if (recordTiming)
            {
                multiplyMilliseconds += ElapsedMilliseconds(start);
            }
#endif
            return;
        }

        output.resize(dimension);

        const float* inputValues = input.data();
        float* outputValues = output.data();
        const int* rowStarts = rowStart.data();
        const int* columnIndices = colIndex.data();
        const Mat3* blockValues = values.data();

        for (int rowBlock = 0; rowBlock < blockCount; ++rowBlock)
        {
            float rowX = 0.0f;
            float rowY = 0.0f;
            float rowZ = 0.0f;
            const int rowBegin = rowStarts[rowBlock];
            const int rowEnd = rowStarts[rowBlock + 1];
            for (int blockIndex = rowBegin; blockIndex < rowEnd; ++blockIndex)
            {
                const int columnBase = columnIndices[blockIndex] * 3;
                const float x = inputValues[columnBase + 0];
                const float y = inputValues[columnBase + 1];
                const float z = inputValues[columnBase + 2];
                const Mat3& block = blockValues[blockIndex];
                const Vec3& column0 = block.columns[0];
                const Vec3& column1 = block.columns[1];
                const Vec3& column2 = block.columns[2];

                rowX += column0.x * x + column1.x * y + column2.x * z;
                rowY += column0.y * x + column1.y * y + column2.y * z;
                rowZ += column0.z * x + column1.z * y + column2.z * z;
            }

            const int rowBase = rowBlock * 3;
            outputValues[rowBase + 0] = rowX;
            outputValues[rowBase + 1] = rowY;
            outputValues[rowBase + 2] = rowZ;
        }

#if defined(PHYSIK_ENABLE_PERF_LOGGING)
        if (recordTiming)
        {
            multiplyMilliseconds += ElapsedMilliseconds(start);
        }
#endif
    }

    int SparseBlockMatrix::FindBlockIndex(int rowBlock, int colBlock) const
    {
        const auto it = blockLookup.find(MakeKey(rowBlock, colBlock));
        if (it == blockLookup.end())
        {
            return -1;
        }

        return it->second;
    }

    std::uint64_t SparseBlockMatrix::MakeKey(int rowBlock, int colBlock)
    {
        return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(rowBlock)) << 32u) |
            static_cast<std::uint32_t>(colBlock);
    }
}
