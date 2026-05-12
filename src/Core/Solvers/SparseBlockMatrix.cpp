#include "PhysiK/Core/Solvers/SparseBlockMatrix.h"

#include <algorithm>
#include <cmath>

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

        Mat3 MassBlock(float mass)
        {
            return Mat3::FromColumns(
                Vec3{mass, 0.0f, 0.0f},
                Vec3{0.0f, mass, 0.0f},
                Vec3{0.0f, 0.0f, mass});
        }

        bool IsValidCoordinate(int rowNode, int colNode, int nodeCount)
        {
            return rowNode >= 0 && rowNode < nodeCount &&
                colNode >= 0 && colNode < nodeCount;
        }
    }

    void SparseBlockMatrix::Clear()
    {
        nodeCount = 0;
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
        int newNodeCount,
        const std::vector<std::pair<int, int>>& blockCoordinates)
    {
        Clear();
        if (newNodeCount <= 0)
        {
            return;
        }

        nodeCount = newNodeCount;
        std::vector<std::vector<int>> rows(static_cast<std::size_t>(nodeCount));
        for (const std::pair<int, int>& coordinate : blockCoordinates)
        {
            const int rowNode = coordinate.first;
            const int colNode = coordinate.second;
            if (!IsValidCoordinate(rowNode, colNode, nodeCount))
            {
                continue;
            }

            rows[static_cast<std::size_t>(rowNode)].push_back(colNode);
        }

        rowStart.assign(static_cast<std::size_t>(nodeCount + 1), 0);
        for (int row = 0; row < nodeCount; ++row)
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

    void SparseBlockMatrix::BuildFromTetConnectivity(
        int newNodeCount,
        const std::vector<Tet>& tets)
    {
        std::vector<std::pair<int, int>> blockCoordinates;
        blockCoordinates.reserve(tets.size() * 16u);

        for (const Tet& tet : tets)
        {
            const int nodes[4] = {tet.node0, tet.node1, tet.node2, tet.node3};
            for (int row = 0; row < 4; ++row)
            {
                for (int column = 0; column < 4; ++column)
                {
                    blockCoordinates.push_back({nodes[row], nodes[column]});
                }
            }
        }

        BuildPattern(newNodeCount, blockCoordinates);
    }

    bool SparseBlockMatrix::AddBlock(int rowNode, int colNode, const Mat3& block)
    {
        const int blockIndex = FindBlockIndex(rowNode, colNode);
        if (blockIndex < 0)
        {
            return false;
        }

        values[static_cast<std::size_t>(blockIndex)] =
            Add(values[static_cast<std::size_t>(blockIndex)], block);
        return true;
    }

    bool SparseBlockMatrix::AddMassToDiagonal(int node, float mass)
    {
        if (!std::isfinite(mass) || mass <= 0.0f)
        {
            return false;
        }

        return AddBlock(node, node, MassBlock(mass));
    }

    void SparseBlockMatrix::AddMassToDiagonal(const std::vector<float>& masses)
    {
        const int count = std::min(nodeCount, static_cast<int>(masses.size()));
        for (int node = 0; node < count; ++node)
        {
            AddMassToDiagonal(node, masses[static_cast<std::size_t>(node)]);
        }
    }

    void SparseBlockMatrix::Multiply(
        const std::vector<float>& input,
        std::vector<float>& output) const
    {
        const std::size_t dimension = static_cast<std::size_t>(std::max(0, nodeCount) * 3);
        output.assign(dimension, 0.0f);
        if (input.size() < dimension || rowStart.size() != static_cast<std::size_t>(nodeCount + 1))
        {
            return;
        }

        for (int rowNode = 0; rowNode < nodeCount; ++rowNode)
        {
            Vec3 rowValue;
            const int rowBegin = rowStart[static_cast<std::size_t>(rowNode)];
            const int rowEnd = rowStart[static_cast<std::size_t>(rowNode + 1)];
            for (int blockIndex = rowBegin; blockIndex < rowEnd; ++blockIndex)
            {
                const int columnNode = colIndex[static_cast<std::size_t>(blockIndex)];
                const std::size_t columnBase = static_cast<std::size_t>(columnNode * 3);
                const Vec3 columnValue{
                    input[columnBase + 0],
                    input[columnBase + 1],
                    input[columnBase + 2]};
                rowValue += values[static_cast<std::size_t>(blockIndex)] * columnValue;
            }

            const std::size_t rowBase = static_cast<std::size_t>(rowNode * 3);
            output[rowBase + 0] = rowValue.x;
            output[rowBase + 1] = rowValue.y;
            output[rowBase + 2] = rowValue.z;
        }
    }

    int SparseBlockMatrix::FindBlockIndex(int rowNode, int colNode) const
    {
        const auto it = blockLookup.find(MakeKey(rowNode, colNode));
        if (it == blockLookup.end())
        {
            return -1;
        }

        return it->second;
    }

    std::uint64_t SparseBlockMatrix::MakeKey(int rowNode, int colNode)
    {
        return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(rowNode)) << 32u) |
            static_cast<std::uint32_t>(colNode);
    }
}
