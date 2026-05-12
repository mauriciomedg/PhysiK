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
    };

    struct ConjugateGradientResult
    {
        int iterations = 0;
        float residualNorm = 0.0f;
        bool converged = false;
    };

    PHYSIK_API ConjugateGradientResult SolveConjugateGradient(
        const SparseBlockMatrix& matrix,
        const std::vector<float>& rhs,
        std::vector<float>& solution,
        const ConjugateGradientSettings& settings = ConjugateGradientSettings{});
}
