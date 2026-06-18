#pragma once

#include <vector>

#include "PhysiK/API/PhysiKAPI.h"
#include "PhysiK/Core/Solvers/Linear/ConjugateGradientSolver.h"
#include "PhysiK/Math/Mat3.h"
#include "PhysiK/Math/SparseBlockMatrix.h"
#include "PhysiK/Math/Vec3.h"

namespace PhysiK
{
    struct LinearSolveSettings
    {
        int maxIterations = 128;
        float tolerance = 1.0e-5f;
        bool useJacobiPreconditioner = true;
    };

    struct LinearSolveResult
    {
        int iterations = 0;
        float residualNorm = 0.0f;
        bool converged = false;
    };

#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
    struct LinearSolverProfileData
    {
        double totalSolveMilliseconds = 0.0;
        double sparseMatrixMultiplyMilliseconds = 0.0;
        double dotProductMilliseconds = 0.0;
        double vectorUpdateMilliseconds = 0.0;
        double preconditionerSetupMilliseconds = 0.0;
        double preconditionerApplyMilliseconds = 0.0;
        int iterations = 0;
        float residualNorm = 0.0f;
        bool converged = false;
    };

    PHYSIK_API void ResetCurrentLinearSolverProfile();
    PHYSIK_API LinearSolverProfileData GetCurrentLinearSolverProfile();
#endif

    class PHYSIK_API LinearSolver
    {
    public:
        virtual ~LinearSolver() = default;

        virtual LinearSolveResult Solve(
            const SparseBlockMatrix& matrix,
            const std::vector<Vec3>& rhs,
            std::vector<Vec3>& solution,
            const LinearSolveSettings& settings) = 0;
    };

    class PHYSIK_API CurrentLinearSolver final : public LinearSolver
    {
    public:
        LinearSolveResult Solve(
            const SparseBlockMatrix& matrix,
            const std::vector<Vec3>& rhs,
            std::vector<Vec3>& solution,
            const LinearSolveSettings& settings) override;

    private:
        std::vector<Mat3> inversePreconditioner;
        std::vector<Vec3> residual;
        std::vector<Vec3> direction;
        std::vector<Vec3> temp;
    };

    PHYSIK_API LinearSolver& GetCurrentLinearSolver();
}
