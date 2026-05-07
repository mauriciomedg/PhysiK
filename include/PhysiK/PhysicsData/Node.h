#pragma once

#include "PhysiK/Math/Vec3.h"

namespace PhysiK
{
    struct Node
    {
        Vec3 position;
        Vec3 velocity;
        Vec3 force;
        float inverseMass = 1.0f;
    };
}
