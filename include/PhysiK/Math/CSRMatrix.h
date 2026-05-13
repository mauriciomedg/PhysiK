#pragma once

#include <vector>

#include "PhysiK/API/PhysiKAPI.h"

namespace PhysiK
{
    struct PHYSIK_API CSRMatrix
    {
        int rowCount = 0;
        int colCount = 0;
        std::vector<int> rowOffsets;
        std::vector<int> columnIndices;
        std::vector<double> values;

        void Clear();
        void Reserve(int nonZeroCount);
        bool IsValid() const;
        void Multiply(
            const std::vector<double>& input,
            std::vector<double>& output) const;
    };
}
