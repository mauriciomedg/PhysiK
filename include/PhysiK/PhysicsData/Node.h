#pragma once

#include "PhysiK/Math/Vec3.h"

namespace PhysiK
{
    struct Node
    {
        Vec3 position;
        Vec3 restPosition;
        Vec3 velocity;
        bool fixed = false;
    };
}
