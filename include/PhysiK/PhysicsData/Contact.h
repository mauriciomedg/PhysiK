#pragma once

#include "PhysiK/API/Handles.h"
#include "PhysiK/Math/Vec3.h"
#include "PhysiK/Math/Vec4.h"

namespace PhysiK
{
    struct Contact
    {
        ComponentHandle sourceComponent;
        ComponentHandle targetComponent;

        int node0 = -1;
        int node1 = -1;
        int node2 = -1;
        int node3 = -1;

        Vec4 barycentric;

        Vec3 worldPoint;
        Vec3 normal;

        float penetrationDepth = 0.0f;
        float stiffness = 0.0f;
        float damping = 0.0f;
    };
}
