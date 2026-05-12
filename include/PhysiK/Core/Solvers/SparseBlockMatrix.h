#pragma once

#include <cstdint>
#include <unordered_map>
#include <utility>
#include <vector>

#include "PhysiK/API/PhysiKAPI.h"
#include "PhysiK/Math/Mat3.h"
#include "PhysiK/PhysicsData/Tet.h"

namespace PhysiK
{
    class PHYSIK_API SparseBlockMatrix
    {
    public:
        int nodeCount = 0;
        std::vector<int> rowStart;
        std::vector<int> colIndex;
        std::vector<Mat3> values;

        void Clear();
        void ClearValues();
        void BuildPattern(
            int newNodeCount,
            const std::vector<std::pair<int, int>>& blockCoordinates);
        void BuildFromTetConnectivity(
            int newNodeCount,
            const std::vector<Tet>& tets);

        bool AddBlock(int rowNode, int colNode, const Mat3& block);
        bool AddMassToDiagonal(int node, float mass);
        void AddMassToDiagonal(const std::vector<float>& masses);
        void Multiply(
            const std::vector<float>& input,
            std::vector<float>& output) const;

        int FindBlockIndex(int rowNode, int colNode) const;

    private:
        std::unordered_map<std::uint64_t, int> blockLookup;

        static std::uint64_t MakeKey(int rowNode, int colNode);
    };
}
