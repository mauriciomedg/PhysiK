#pragma once

#if defined(PHYSIK_ENABLE_MKL)

#include <vector>

#include "PhysiK/API/PhysiKAPI.h"
#include "PhysiK/Core/Solvers/Linear/LinearSolver.h"

namespace PhysiK
{
    class PHYSIK_API MKLLinearSolver final : public LinearSolver
    {
    public:
        LinearSolveResult Solve(
            const SparseBlockMatrix& matrix,
            const std::vector<float>& rhs,
            std::vector<float>& solution,
            const LinearSolveSettings& settings) override;
    };
}

#endif
