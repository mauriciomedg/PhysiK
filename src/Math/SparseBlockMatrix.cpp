#include "PhysiK/Math/SparseBlockMatrix.h"

#include <algorithm>

namespace PhysiK
{
    namespace
    {
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
        const std::size_t dimension = static_cast<std::size_t>(std::max(0, blockCount) * 3);
        output.assign(dimension, 0.0f);
        if (input.size() < dimension || rowStart.size() != static_cast<std::size_t>(blockCount + 1))
        {
            return;
        }

        for (int rowBlock = 0; rowBlock < blockCount; ++rowBlock)
        {
            Vec3 rowValue;
            const int rowBegin = rowStart[static_cast<std::size_t>(rowBlock)];
            const int rowEnd = rowStart[static_cast<std::size_t>(rowBlock + 1)];
            for (int blockIndex = rowBegin; blockIndex < rowEnd; ++blockIndex)
            {
                const int columnBlock = colIndex[static_cast<std::size_t>(blockIndex)];
                const std::size_t columnBase = static_cast<std::size_t>(columnBlock * 3);
                const Vec3 columnValue{
                    input[columnBase + 0],
                    input[columnBase + 1],
                    input[columnBase + 2]};
                rowValue += values[static_cast<std::size_t>(blockIndex)] * columnValue;
            }

            const std::size_t rowBase = static_cast<std::size_t>(rowBlock * 3);
            output[rowBase + 0] = rowValue.x;
            output[rowBase + 1] = rowValue.y;
            output[rowBase + 2] = rowValue.z;
        }
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
