#pragma once

#include "PhysiK/Math/Mat3.h"

namespace PhysiK
{
    struct Tet
    {
        int node0 = -1;
        int node1 = -1;
        int node2 = -1;
        int node3 = -1;

        float restVolume = 0.0f;
        Mat3 restDmInverse;
        float stiffness = 25.0f;
        float damping = 0.25f;
    };
}
