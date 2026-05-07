#pragma once

#include <algorithm>

#include "PhysiK/Math/Vec3.h"

namespace PhysiK
{
    struct AABB
    {
        Vec3 minimum;
        Vec3 maximum;

        bool Contains(const Vec3& point) const
        {
            return point.x >= minimum.x && point.x <= maximum.x &&
                point.y >= minimum.y && point.y <= maximum.y &&
                point.z >= minimum.z && point.z <= maximum.z;
        }

        void ExpandToInclude(const Vec3& point)
        {
            minimum.x = std::min(minimum.x, point.x);
            minimum.y = std::min(minimum.y, point.y);
            minimum.z = std::min(minimum.z, point.z);
            maximum.x = std::max(maximum.x, point.x);
            maximum.y = std::max(maximum.y, point.y);
            maximum.z = std::max(maximum.z, point.z);
        }
    };
}
