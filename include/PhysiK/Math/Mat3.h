#pragma once

#include <cmath>

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

        static Mat3 Zero()
        {
            return FromColumns(Vec3{}, Vec3{}, Vec3{});
        }

        static Mat3 FromColumns(const Vec3& column0, const Vec3& column1, const Vec3& column2)
        {
            Mat3 matrix;
            matrix.columns[0] = column0;
            matrix.columns[1] = column1;
            matrix.columns[2] = column2;
            return matrix;
        }

        Vec3 operator*(const Vec3& vector) const
        {
            return columns[0] * vector.x + columns[1] * vector.y + columns[2] * vector.z;
        }

        Mat3 operator-(const Mat3& other) const
        {
            return FromColumns(
                columns[0] - other.columns[0],
                columns[1] - other.columns[1],
                columns[2] - other.columns[2]);
        }
    };

    inline Mat3 operator*(const Mat3& a, const Mat3& b)
    {
        return Mat3::FromColumns(a * b.columns[0], a * b.columns[1], a * b.columns[2]);
    }

    inline Mat3 Transpose(const Mat3& matrix)
    {
        return Mat3::FromColumns(
            Vec3{matrix.columns[0].x, matrix.columns[1].x, matrix.columns[2].x},
            Vec3{matrix.columns[0].y, matrix.columns[1].y, matrix.columns[2].y},
            Vec3{matrix.columns[0].z, matrix.columns[1].z, matrix.columns[2].z});
    }

    inline float Determinant(const Mat3& matrix)
    {
        return Dot(matrix.columns[0], Cross(matrix.columns[1], matrix.columns[2]));
    }

    inline Mat3 Inverse(const Mat3& matrix)
    {
        const Vec3 row0 = Cross(matrix.columns[1], matrix.columns[2]);
        const Vec3 row1 = Cross(matrix.columns[2], matrix.columns[0]);
        const Vec3 row2 = Cross(matrix.columns[0], matrix.columns[1]);
        const float determinant = Dot(matrix.columns[0], row0);

        if (std::abs(determinant) <= 0.000001f)
        {
            return Mat3{};
        }

        const float inverseDeterminant = 1.0f / determinant;
        return Transpose(Mat3::FromColumns(
            row0 * inverseDeterminant,
            row1 * inverseDeterminant,
            row2 * inverseDeterminant));
    }
}
