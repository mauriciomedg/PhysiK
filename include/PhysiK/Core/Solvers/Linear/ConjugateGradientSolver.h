#pragma once

#include <vector>

#include "PhysiK/API/PhysiKAPI.h"
#include "PhysiK/Math/Mat3.h"
#include "PhysiK/Math/SparseBlockMatrix.h"
#include "PhysiK/Math/Vec3.h"

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

    PHYSIK_API ConjugateGradientResult SolvePreconditionedConjugateGradient(
        std::vector<Vec3>& x,
        const SparseBlockMatrix& A,
        const std::vector<Vec3>& b,
        int maxIterations,
        float tolerance,
        const std::vector<Mat3>& MInv,
        std::vector<Vec3>& r,
        std::vector<Vec3>& d,
        std::vector<Vec3>& qOrS);
}
