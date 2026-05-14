#pragma once

#include <vector>

#include "PhysiK/API/PhysiKAPI.h"
#include "PhysiK/Math/SparseBlockMatrix.h"

namespace PhysiK
{
    enum class LinearSolverBackend
    {
        Current,
        MKL
    };

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

    class PHYSIK_API LinearSolver
    {
    public:
        virtual ~LinearSolver() = default;

        virtual LinearSolveResult Solve(
            const SparseBlockMatrix& matrix,
            const std::vector<float>& rhs,
            std::vector<float>& solution,
            const LinearSolveSettings& settings) = 0;
    };

    class PHYSIK_API CurrentLinearSolver final : public LinearSolver
    {
    public:
        LinearSolveResult Solve(
            const SparseBlockMatrix& matrix,
            const std::vector<float>& rhs,
            std::vector<float>& solution,
            const LinearSolveSettings& settings) override;
    };

    PHYSIK_API bool IsLinearSolverBackendAvailable(LinearSolverBackend backend);

    PHYSIK_API LinearSolver& GetLinearSolver(LinearSolverBackend backend);
}
