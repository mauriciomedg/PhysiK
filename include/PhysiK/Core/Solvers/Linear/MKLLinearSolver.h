#pragma once

#if defined(PHYSIK_ENABLE_MKL)

#include <vector>

#include "PhysiK/API/PhysiKAPI.h"
#include "PhysiK/Core/Solvers/Linear/LinearSolver.h"
#include "PhysiK/Math/CSRMatrix.h"

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

        LinearSolveResult SolveSPD(
            const CSRMatrix& matrix,
            const std::vector<double>& rhs,
            std::vector<double>& solution,
            const LinearSolveSettings& settings);
    };
}

#endif
