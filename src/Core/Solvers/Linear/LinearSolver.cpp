#include "PhysiK/Core/Solvers/Linear/LinearSolver.h"

#include <algorithm>
#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
#include <chrono>
#endif

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

#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
        using Clock = std::chrono::steady_clock;

        LinearSolverProfileData currentProfile;

        double MillisecondsBetween(Clock::time_point start, Clock::time_point end)
        {
            return std::chrono::duration<double, std::milli>(end - start).count();
        }
#endif
    }

#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
    void ResetCurrentLinearSolverProfile()
    {
        currentProfile = LinearSolverProfileData{};
        ResetConjugateGradientProfile();
    }

    LinearSolverProfileData GetCurrentLinearSolverProfile()
    {
        return currentProfile;
    }
#endif

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

#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
        ResetCurrentLinearSolverProfile();
        const Clock::time_point totalStart = Clock::now();
#endif
        const ConjugateGradientResult cgResult =
            SolveConjugateGradient(matrix, rhs, solution, cgSettings);
#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
        const double totalMilliseconds = MillisecondsBetween(totalStart, Clock::now());
        const ConjugateGradientProfileData cgProfile = GetConjugateGradientProfile();
        currentProfile.totalSolveMilliseconds = totalMilliseconds;
        currentProfile.sparseMatrixMultiplyMilliseconds =
            cgProfile.sparseMatrixMultiplyMilliseconds;
        currentProfile.dotProductMilliseconds = cgProfile.dotProductMilliseconds;
        currentProfile.vectorUpdateMilliseconds = cgProfile.vectorUpdateMilliseconds;
        currentProfile.preconditionerSetupMilliseconds =
            cgProfile.preconditionerSetupMilliseconds;
        currentProfile.preconditionerApplyMilliseconds =
            cgProfile.preconditionerApplyMilliseconds;
        currentProfile.iterations = cgResult.iterations;
        currentProfile.residualNorm = cgResult.residualNorm;
        currentProfile.converged = cgResult.converged;
#endif

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
