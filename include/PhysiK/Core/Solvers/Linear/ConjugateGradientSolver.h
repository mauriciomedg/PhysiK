#pragma once

#include <vector>

#include "PhysiK/API/PhysiKAPI.h"
#include "PhysiK/Math/SparseBlockMatrix.h"

namespace PhysiK
{
    struct ConjugateGradientSettings
    {
        int maxIterations = 128;
        float tolerance = 1.0e-5f;
        bool useJacobiPreconditioner = true;
        SparseBlockMatrixMultiplyMode multiplyMode = SparseBlockMatrixMultiplyMode::Serial;
        int multiplyWorkerCount = GetSparseBlockMatrixMultiplyWorkerCount();
    };

    struct ConjugateGradientResult
    {
        int iterations = 0;
        float residualNorm = 0.0f;
        bool converged = false;
    };

    struct ConjugateGradientScratch
    {
        std::vector<float> residual;
        std::vector<float> direction;
        std::vector<float> matrixDirection;
        std::vector<float> preconditionedResidual;
        std::vector<float> inverseDiagonal;
        std::vector<int> diagonalBlockIndices;
        std::vector<int> cachedRowStart;
        std::vector<int> cachedColumnIndex;
        int cachedBlockCount = -1;
    };

#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
    struct ConjugateGradientProfileData
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

    PHYSIK_API void ResetConjugateGradientProfile();
    PHYSIK_API ConjugateGradientProfileData GetConjugateGradientProfile();
#endif

    PHYSIK_API ConjugateGradientResult SolveConjugateGradient(
        const SparseBlockMatrix& matrix,
        const std::vector<float>& rhs,
        std::vector<float>& solution,
        const ConjugateGradientSettings& settings = ConjugateGradientSettings{});

    PHYSIK_API ConjugateGradientResult SolveConjugateGradient(
        const SparseBlockMatrix& matrix,
        const std::vector<float>& rhs,
        std::vector<float>& solution,
        ConjugateGradientScratch& scratch,
        const ConjugateGradientSettings& settings = ConjugateGradientSettings{});
}
