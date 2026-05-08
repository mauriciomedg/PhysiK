#pragma once

#include <vector>

namespace PhysiK
{
    class CollisionDetectionEngine;
    struct Contact;
    class SolverData;
    class World;

    class Component
    {
    public:
        bool active = true;

        virtual void UpdateFrame(World& world, float dt)
        {
            (void)world;
            (void)dt;
        }

        virtual void QueryContacts(
            World& world,
            CollisionDetectionEngine& collisionDetectionEngine,
            std::vector<Contact>& outContacts)
        {
            (void)world;
            (void)collisionDetectionEngine;
            (void)outContacts;
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
