#pragma once

namespace PhysiK
{
    class SolverData;
    class World;

    class PhysicsConnection
    {
    public:
        virtual void UpdateSystem(World& world, SolverData& solverData, float dt) = 0;
        virtual ~PhysicsConnection() = default;
    };
}
