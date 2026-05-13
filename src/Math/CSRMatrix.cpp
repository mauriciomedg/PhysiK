#include "PhysiK/Math/CSRMatrix.h"

#include <algorithm>

namespace PhysiK
{
    void CSRMatrix::Clear()
    {
        rowCount = 0;
        colCount = 0;
        rowOffsets.clear();
        columnIndices.clear();
        values.clear();
    }

    void CSRMatrix::Reserve(int nonZeroCount)
    {
        const std::size_t capacity =
            static_cast<std::size_t>(std::max(0, nonZeroCount));
        columnIndices.reserve(capacity);
        values.reserve(capacity);
    }

    bool CSRMatrix::IsValid() const
    {
        if (rowCount < 0 || colCount < 0)
        {
            return false;
        }

        if (rowOffsets.size() != static_cast<std::size_t>(rowCount + 1))
        {
            return false;
        }

        if (rowOffsets.empty() || rowOffsets.front() != 0)
        {
            return false;
        }

        if (columnIndices.size() != values.size())
        {
            return false;
        }

        for (int row = 0; row < rowCount; ++row)
        {
            const int rowBegin = rowOffsets[static_cast<std::size_t>(row)];
            const int rowEnd = rowOffsets[static_cast<std::size_t>(row + 1)];
            if (rowBegin > rowEnd)
            {
                return false;
            }
        }

        if (rowOffsets.back() != static_cast<int>(values.size()))
        {
            return false;
        }

        for (int column : columnIndices)
        {
            if (column < 0 || column >= colCount)
            {
                return false;
            }
        }

        return true;
    }

    void CSRMatrix::Multiply(
        const std::vector<double>& input,
        std::vector<double>& output) const
    {
        output.assign(static_cast<std::size_t>(std::max(0, rowCount)), 0.0);
        if (!IsValid() || input.size() < static_cast<std::size_t>(std::max(0, colCount)))
        {
            return;
        }

        for (int row = 0; row < rowCount; ++row)
        {
            double sum = 0.0;
            const int rowBegin = rowOffsets[static_cast<std::size_t>(row)];
            const int rowEnd = rowOffsets[static_cast<std::size_t>(row + 1)];
            for (int valueIndex = rowBegin; valueIndex < rowEnd; ++valueIndex)
            {
                const std::size_t index = static_cast<std::size_t>(valueIndex);
                const int column = columnIndices[index];
                sum += values[index] * input[static_cast<std::size_t>(column)];
            }

            output[static_cast<std::size_t>(row)] = sum;
        }
    }
}
