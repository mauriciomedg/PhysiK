#pragma once

#include "PhysiK/Math/Vec3.h"
#include "PhysiK/Math/Vec4.h"

namespace PhysiK
{
    struct Contact
    {
        int tetNode0 = -1;
        int tetNode1 = -1;
        int tetNode2 = -1;
        int tetNode3 = -1;

        Vec4 barycentric;

        Vec3 worldPoint;
        Vec3 normal;

        float penetrationDepth = 0.0f;
        float stiffness = 0.0f;
        float damping = 0.0f;
    };
}
