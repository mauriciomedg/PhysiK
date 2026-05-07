#pragma once

#include "PhysiK/Math/Vec3.h"

namespace PhysiK
{
    struct Mat3
    {
        Vec3 columns[3] = {
            Vec3{1.0f, 0.0f, 0.0f},
            Vec3{0.0f, 1.0f, 0.0f},
            Vec3{0.0f, 0.0f, 1.0f}};

        static Mat3 Identity()
        {
            return Mat3{};
        }

        Vec3 operator*(const Vec3& vector) const
        {
            return columns[0] * vector.x + columns[1] * vector.y + columns[2] * vector.z;
        }
    };
}
