#include "PhysiK/Core/Solvers/Linear/LinearSolver.h"

#include <algorithm>

#include "PhysiK/Core/Solvers/Linear/ConjugateGradientSolver.h"

namespace PhysiK
{
    namespace
    {
        CurrentLinearSolver& CurrentSolver()
        {
            static CurrentLinearSolver solver;
            return solver;
        }
    }

    LinearSolveResult CurrentLinearSolver::Solve(
        const SparseBlockMatrix& matrix,
        const std::vector<float>& rhs,
        std::vector<float>& solution,
        const LinearSolveSettings& settings)
    {
        ConjugateGradientSettings cgSettings;
        cgSettings.maxIterations = std::max(1, settings.maxIterations);
        cgSettings.tolerance = settings.tolerance;
        cgSettings.useJacobiPreconditioner = settings.useJacobiPreconditioner;

        const ConjugateGradientResult cgResult =
            SolveConjugateGradient(matrix, rhs, solution, cgSettings);

        LinearSolveResult result;
        result.iterations = cgResult.iterations;
        result.residualNorm = cgResult.residualNorm;
        result.converged = cgResult.converged;
        return result;
    }

    LinearSolver& GetCurrentLinearSolver()
    {
        return CurrentSolver();
    }
}
