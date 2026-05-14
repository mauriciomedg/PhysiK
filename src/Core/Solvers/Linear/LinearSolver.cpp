#include "PhysiK/Core/Solvers/Linear/LinearSolver.h"

#include <algorithm>

#include "PhysiK/Core/Solvers/Linear/ConjugateGradientSolver.h"

#if defined(PHYSIK_ENABLE_MKL)
#include "PhysiK/Core/Solvers/Linear/MKLLinearSolver.h"
#endif

namespace PhysiK
{
    namespace
    {
        CurrentLinearSolver& CurrentSolver()
        {
            static CurrentLinearSolver solver;
            return solver;
        }

#if defined(PHYSIK_ENABLE_MKL)
        MKLLinearSolver& MKLSolver()
        {
            static MKLLinearSolver solver;
            return solver;
        }
#endif
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

    bool IsLinearSolverBackendAvailable(LinearSolverBackend backend)
    {
        switch (backend)
        {
        case LinearSolverBackend::Current:
            return true;
        case LinearSolverBackend::MKL:
#if defined(PHYSIK_ENABLE_MKL)
            return true;
#else
            return false;
#endif
        default:
            return false;
        }
    }

    LinearSolver& GetLinearSolver(LinearSolverBackend backend)
    {
        if (backend == LinearSolverBackend::MKL)
        {
#if defined(PHYSIK_ENABLE_MKL)
            return MKLSolver();
#else
            return CurrentSolver();
#endif
        }

        return CurrentSolver();
    }
}
