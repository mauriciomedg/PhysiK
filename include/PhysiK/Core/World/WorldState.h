#pragma once

#include <vector>

#include "PhysiK/Math/Mat3.h"
#include "PhysiK/Math/Quaternion.h"
#include "PhysiK/Math/Vec3.h"

namespace PhysiK
{
    struct WorldState
    {
        std::vector<Vec3> positions;
        std::vector<Vec3> velocities;
        std::vector<Vec3> forces;
        std::vector<float> masses;
        std::vector<Quaternion> orientations;
        std::vector<Vec3> angularVelocities;
        std::vector<Vec3> torques;
        std::vector<Mat3> inverseInertias;

        void Clear()
        {
            positions.clear();
            velocities.clear();
            forces.clear();
            masses.clear();
            orientations.clear();
            angularVelocities.clear();
            torques.clear();
            inverseInertias.clear();
        }

        int AddNodeState(
            const Vec3& position,
            const Vec3& velocity = Vec3{},
            float mass = 0.0f)
        {
            const int index = static_cast<int>(positions.size());
            positions.push_back(position);
            velocities.push_back(velocity);
            forces.push_back(Vec3{});
            masses.push_back(mass);
            orientations.push_back(Quaternion::Identity());
            angularVelocities.push_back(Vec3{});
            torques.push_back(Vec3{});
            inverseInertias.push_back(Mat3::Zero());
            return index;
        }

        bool IsValidStateIndex(int index) const
        {
            return index >= 0 &&
                index < static_cast<int>(positions.size()) &&
                index < static_cast<int>(velocities.size()) &&
                index < static_cast<int>(forces.size()) &&
                index < static_cast<int>(masses.size()) &&
                index < static_cast<int>(orientations.size()) &&
                index < static_cast<int>(angularVelocities.size()) &&
                index < static_cast<int>(torques.size()) &&
                index < static_cast<int>(inverseInertias.size());
        }
    };
}
