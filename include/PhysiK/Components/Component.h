#pragma once

#include <vector>

#include "PhysiK/Core/Events/PhysicsEvent.h"

namespace PhysiK
{
    class SolverData;
    class World;

    class Component
    {
    public:
        bool active = true;
        std::vector<PhysicsEventType> listenedEvents;
        std::vector<PhysicsEventType> emittedEvents;

        virtual void PreUpdate(World&, float)
        {
        }

        virtual void UpdateSystem(World&, SolverData&, float)
        {
        }

        virtual void PostUpdate(World&, float)
        {
        }

        virtual void OnPhysicsEvent(const PhysicsEvent&)
        {
        }

        virtual ~Component() = default;
    };
}
