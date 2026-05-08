#pragma once

namespace PhysiK
{
    class SolverData;
    class World;

    class PhysicsModel
    {
    public:
        virtual void UpdateSystem(World& world, SolverData& solverData, float dt) = 0;
        virtual ~PhysicsModel() = default;
    };
}
