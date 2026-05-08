#pragma once

namespace PhysiK
{
    class SolverData;
    class World;

    class Component
    {
    public:
        bool active = true;

        virtual void Update(World& world, float dt)
        {
            (void)world;
            (void)dt;
        }

        virtual void UpdateSystem(World& world, SolverData& solverData, float dt)
        {
            (void)world;
            (void)solverData;
            (void)dt;
        }

        virtual ~Component() = default;
    };
}
