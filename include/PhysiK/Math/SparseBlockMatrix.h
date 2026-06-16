#pragma once

#include <cstdint>
#include <unordered_map>
#include <utility>
#include <vector>

#include "PhysiK/API/PhysiKAPI.h"
#include "PhysiK/Math/Mat3.h"

namespace PhysiK
{
    PHYSIK_API void SetSparseBlockMatrixTimingEnabled(bool enabled);
    PHYSIK_API void ResetSparseBlockMatrixTiming();
    PHYSIK_API double GetSparseBlockMatrixMultiplyMilliseconds();

    class PHYSIK_API SparseBlockMatrix
    {
    public:
        int blockCount = 0;
        std::vector<int> rowStart;
        std::vector<int> colIndex;
        std::vector<Mat3> values;

        void Clear();
        void ClearValues();
        void BuildPattern(
            int newBlockCount,
            const std::vector<std::pair<int, int>>& blockCoordinates);

        bool AddBlock(int rowBlock, int colBlock, const Mat3& block);
        bool AddBlockAtIndex(int blockIndex, const Mat3& block);
        void Multiply(
            const std::vector<float>& input,
            std::vector<float>& output) const;

        int FindBlockIndex(int rowBlock, int colBlock) const;

    private:
        std::unordered_map<std::uint64_t, int> blockLookup;

        static std::uint64_t MakeKey(int rowBlock, int colBlock);
    };
}
