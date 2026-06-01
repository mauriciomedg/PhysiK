#pragma once

#include <vector>

#include "PhysiK/API/PhysiKAPI.h"
#include "PhysiK/Math/SparseBlockMatrix.h"

namespace PhysiK
{
    struct ConjugateGradientSettings
    {
        int maxIterations = 128;
        float tolerance = 1.0e-4f;
        bool useJacobiPreconditioner = true;
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
