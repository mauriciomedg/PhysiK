#pragma once

#include "PhysiK/Math/Vec3.h"

namespace PhysiK
{
    struct Quaternion
    {
        float w = 1.0f;
        float x = 0.0f;
        float y = 0.0f;
        float z = 0.0f;

        constexpr Quaternion() = default;
        constexpr Quaternion(float wValue, float xValue, float yValue, float zValue)
            : w(wValue), x(xValue), y(yValue), z(zValue)
        {
        }

        static Quaternion Identity()
        {
            return Quaternion{};
        }

        Quaternion Conjugate() const
        {
            return Quaternion{w, -x, -y, -z};
        }

        Quaternion operator*(const Quaternion& other) const
        {
            return Quaternion{
                w * other.w - x * other.x - y * other.y - z * other.z,
                w * other.x + x * other.w + y * other.z - z * other.y,
                w * other.y - x * other.z + y * other.w + z * other.x,
                w * other.z + x * other.y - y * other.x + z * other.w};
        }

        Vec3 operator*(const Vec3& vector) const
        {
            const Quaternion rotated = (*this) * Quaternion{0.0f, vector.x, vector.y, vector.z} * Conjugate();
            return Vec3{rotated.x, rotated.y, rotated.z};
        }
    };
}
