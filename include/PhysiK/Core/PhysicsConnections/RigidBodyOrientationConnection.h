#pragma once

#include "PhysiK/API/Handles.h"
#include "PhysiK/Core/PhysicsConnections/PhysicsConnection.h"
#include "PhysiK/Math/Quaternion.h"

namespace PhysiK
{
    struct RigidBodyOrientationConnection : public PhysicsConnection
    {
        RigidBodyHandle rigidBody;

        Quaternion targetOrientation;

        float stiffness = 0.0f;
        float damping = 0.0f;

        void UpdateSystem(World& world, SolverData& solverData, float dt) override;
    };
}
