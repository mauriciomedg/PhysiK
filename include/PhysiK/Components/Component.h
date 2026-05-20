#pragma once

#include <vector>

#include "PhysiK/Core/Events/PhysicsEvent.h"

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
        std::vector<PhysicsEventType> listenedEvents;
        std::vector<PhysicsEventType> emittedEvents;

        virtual void UpdateFrame(World& world, float dt)
        {
            (void)world;
            (void)dt;
        }

        virtual void UpdateKinematicTarget(World& world)
        {
            (void)world;
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

        virtual void OnPhysicsEvent(const PhysicsEvent& event)
        {
            (void)event;
        }

        virtual ~Component() = default;
    };
}
