#pragma once

#include "PhysiK/API/Handles.h"
#include "PhysiK/Core/PhysicsConnections/PhysicsConnection.h"
#include "PhysiK/Math/Vec3.h"

namespace PhysiK
{
    struct RigidBodyConnection : public PhysicsConnection
    {
        RigidBodyHandle rigidBody;

        Vec3 localPoint;
        Vec3 targetPosition;

        float stiffness = 0.0f;
        float damping = 0.0f;

        void UpdateSystem(World& world, SolverData& solverData, float dt) override;
    };
}
