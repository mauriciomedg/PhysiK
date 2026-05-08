#include "PhysiK/Core/PhysicsConnections/RigidBodyConnection.h"

namespace PhysiK
{
    void RigidBodyConnection::UpdateSystem(World& world, SolverData& solverData, float dt)
    {
        (void)world;
        (void)solverData;
        (void)dt;
        // TODO: Add linear/angular force contributions to SolverData.
        // This connection must not directly modify rigid-body transforms or velocities.
    }
}
