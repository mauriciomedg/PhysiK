#include "PhysiK/Core/Solvers/Linear/MKLLinearSolver.h"

#if defined(PHYSIK_ENABLE_MKL)

namespace PhysiK
{
    LinearSolveResult MKLLinearSolver::Solve(
        const SparseBlockMatrix&,
        const std::vector<float>&,
        std::vector<float>& solution,
        const LinearSolveSettings&)
    {
        solution.clear();
        return LinearSolveResult{};
    }
}

#endif
