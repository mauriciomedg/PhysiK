#pragma once

namespace PhysiK
{
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

        virtual ~Component() = default;
    };
}
