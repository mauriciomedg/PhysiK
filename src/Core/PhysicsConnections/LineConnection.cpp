#include "PhysiK/Core/PhysicsConnections/LineConnection.h"

namespace PhysiK
{
    void LineConnection::UpdateSystem(World& world, SolverData& solverData, float dt)
    {
        (void)world;
        (void)solverData;
        (void)dt;
        // TODO: Add force, damping, and stiffness contributions to SolverData.
        // This connection must not directly modify node positions or velocities.
    }
}
