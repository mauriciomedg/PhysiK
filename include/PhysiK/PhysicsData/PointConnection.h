#pragma once

#include "PhysiK/Math/Vec3.h"
#include "PhysiK/Math/Vec4.h"

namespace PhysiK
{
    struct PointConnection
    {
        int node0 = -1;
        int node1 = -1;
        int node2 = -1;
        int node3 = -1;

        Vec4 barycentric;
        Vec3 targetPosition;

        float stiffness = 0.0f;
        float damping = 0.0f;
    };
}
