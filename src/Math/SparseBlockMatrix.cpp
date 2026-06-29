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

        std::pair<int, int> CanonicalCoordinate(int rowBlock, int colBlock)
        {
            return rowBlock <= colBlock
                ? std::pair<int, int>{rowBlock, colBlock}
                : std::pair<int, int>{colBlock, rowBlock};
        }

        void AddMatVec(
            const Mat3& block,
            const Vec3& input,
            float& rowX,
            float& rowY,
            float& rowZ)
        {
            const float x = input.x;
            const float y = input.y;
            const float z = input.z;
            const Vec3& column0 = block.columns[0];
            const Vec3& column1 = block.columns[1];
            const Vec3& column2 = block.columns[2];

            rowX += column0.x * x + column1.x * y + column2.x * z;
            rowY += column0.y * x + column1.y * y + column2.y * z;
            rowZ += column0.z * x + column1.z * y + column2.z * z;
        }

        void AddTransposeMatVec(
            const Mat3& block,
            const Vec3& input,
            float& rowX,
            float& rowY,
            float& rowZ)
        {
            const float x = input.x;
            const float y = input.y;
            const float z = input.z;
            const Vec3& column0 = block.columns[0];
            const Vec3& column1 = block.columns[1];
            const Vec3& column2 = block.columns[2];

            rowX += column0.x * x + column0.y * y + column0.z * z;
            rowY += column1.x * x + column1.y * y + column1.z * z;
            rowZ += column2.x * x + column2.y * y + column2.z * z;
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
            int rowBlock = coordinate.first;
            int colBlock = coordinate.second;
            if (!IsValidCoordinate(rowBlock, colBlock, blockCount))
            {
                continue;
            }

            const std::pair<int, int> canonicalCoordinate =
                CanonicalCoordinate(rowBlock, colBlock);
            rowBlock = canonicalCoordinate.first;
            colBlock = canonicalCoordinate.second;
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

        const Mat3 contribution =
            rowBlock <= colBlock ? block : Transpose(block);
        values[static_cast<std::size_t>(blockIndex)] =
            Add(values[static_cast<std::size_t>(blockIndex)], contribution);
        return true;
    }

    void SparseBlockMatrix::Multiply(
        const std::vector<float>& input,
        std::vector<float>& output) const
    {
        const std::size_t dimension = static_cast<std::size_t>(std::max(0, blockCount) * 3);
        if (input.size() < dimension || rowStart.size() != static_cast<std::size_t>(blockCount + 1))
        {
            output.assign(dimension, 0.0f);
            return;
        }

        output.resize(dimension);
        std::fill(output.begin(), output.end(), 0.0f);

        const float* inputValues = input.data();
        float* outputValues = output.data();
        const int* rowStarts = rowStart.data();
        const int* columnIndices = colIndex.data();
        const Mat3* blockValues = values.data();

        for (int rowBlock = 0; rowBlock < blockCount; ++rowBlock)
        {
            const int rowBegin = rowStarts[rowBlock];
            const int rowEnd = rowStarts[rowBlock + 1];
            for (int blockIndex = rowBegin; blockIndex < rowEnd; ++blockIndex)
            {
                const int columnBlock = columnIndices[blockIndex];
                const int rowBase = rowBlock * 3;
                const int columnBase = columnBlock * 3;
                const Mat3& block = blockValues[blockIndex];
                const Vec3& column0 = block.columns[0];
                const Vec3& column1 = block.columns[1];
                const Vec3& column2 = block.columns[2];
                const float columnX = inputValues[columnBase + 0];
                const float columnY = inputValues[columnBase + 1];
                const float columnZ = inputValues[columnBase + 2];

                outputValues[rowBase + 0] +=
                    column0.x * columnX + column1.x * columnY + column2.x * columnZ;
                outputValues[rowBase + 1] +=
                    column0.y * columnX + column1.y * columnY + column2.y * columnZ;
                outputValues[rowBase + 2] +=
                    column0.z * columnX + column1.z * columnY + column2.z * columnZ;

                if (rowBlock == columnBlock)
                {
                    continue;
                }

                const float rowX = inputValues[rowBase + 0];
                const float rowY = inputValues[rowBase + 1];
                const float rowZ = inputValues[rowBase + 2];
                outputValues[columnBase + 0] +=
                    column0.x * rowX + column0.y * rowY + column0.z * rowZ;
                outputValues[columnBase + 1] +=
                    column1.x * rowX + column1.y * rowY + column1.z * rowZ;
                outputValues[columnBase + 2] +=
                    column2.x * rowX + column2.y * rowY + column2.z * rowZ;
            }
        }

    }

    void SparseBlockMatrix::Multiply(
        const std::vector<Vec3>& input,
        std::vector<Vec3>& output) const
    {
        const std::size_t blockDimension =
            static_cast<std::size_t>(std::max(0, blockCount));

        if (input.size() < blockDimension ||
            rowStart.size() != static_cast<std::size_t>(blockCount + 1))
        {
            output.assign(blockDimension, Vec3{});
            return;
        }

        if (output.size() != blockDimension)
        {
            output.resize(blockDimension);
        }
        std::fill(output.begin(), output.end(), Vec3{});

        const Vec3* inputValues = input.data();
        Vec3* outputValues = output.data();
        const int* rowStarts = rowStart.data();
        const int* columnIndices = colIndex.data();
        const Mat3* blockValues = values.data();

        for (int rowBlock = 0; rowBlock < blockCount; ++rowBlock)
        {
            const int rowBegin = rowStarts[rowBlock];
            const int rowEnd = rowStarts[rowBlock + 1];

            for (int blockIndex = rowBegin; blockIndex < rowEnd; ++blockIndex)
            {
                const int columnBlock = columnIndices[blockIndex];
                const Mat3& block = blockValues[blockIndex];
                float rowX = outputValues[rowBlock].x;
                float rowY = outputValues[rowBlock].y;
                float rowZ = outputValues[rowBlock].z;
                AddMatVec(block, inputValues[columnBlock], rowX, rowY, rowZ);
                outputValues[rowBlock] = Vec3{rowX, rowY, rowZ};

                if (rowBlock == columnBlock)
                {
                    continue;
                }

                float columnX = outputValues[columnBlock].x;
                float columnY = outputValues[columnBlock].y;
                float columnZ = outputValues[columnBlock].z;
                AddTransposeMatVec(
                    block,
                    inputValues[rowBlock],
                    columnX,
                    columnY,
                    columnZ);
                outputValues[columnBlock] = Vec3{columnX, columnY, columnZ};
            }
        }

    }

    int SparseBlockMatrix::FindBlockIndex(int rowBlock, int colBlock) const
    {
        const std::pair<int, int> canonicalCoordinate =
            CanonicalCoordinate(rowBlock, colBlock);
        const auto it = blockLookup.find(
            MakeKey(canonicalCoordinate.first, canonicalCoordinate.second));
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
