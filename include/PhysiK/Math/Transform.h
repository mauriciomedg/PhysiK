#pragma once

#include "PhysiK/Math/Quaternion.h"
#include "PhysiK/Math/Vec3.h"

namespace PhysiK
{
    struct Transform
    {
        Vec3 position;
        Quaternion rotation;

        Vec3 TransformPoint(const Vec3& localPoint) const
        {
            return position + rotation * localPoint;
        }
    };
}
