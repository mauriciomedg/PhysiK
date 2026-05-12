#pragma once

#include "PhysiK/Math/Mat3.h"
#include "PhysiK/Math/Vec3.h"

namespace PhysiK
{
    struct Tet
    {
        int node0 = -1;
        int node1 = -1;
        int node2 = -1;
        int node3 = -1;

        bool active = true;

        float restVolume = 0.0f;
        Mat3 restDmInverse;
        Vec3 restPositions[4];
        Vec3 shapeFunctionGradients[4];
        float youngModulus = 25.0f;
        float poissonRatio = 0.3f;
        float damping = 0.25f;
    };
}
