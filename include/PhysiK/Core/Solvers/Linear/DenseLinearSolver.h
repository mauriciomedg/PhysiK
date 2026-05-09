#pragma once

#include <vector>

namespace PhysiK
{
    class DenseLinearSolver
    {
    public:
        static bool Solve(
            const std::vector<float>& matrix,
            const std::vector<float>& rhs,
            std::vector<float>& solution,
            int dimension);
    };
}
