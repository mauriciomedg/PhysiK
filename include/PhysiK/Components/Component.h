#pragma once

#include <vector>

#include "PhysiK/Components/ComponentExecutionPriority.h"
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

        virtual ComponentExecutionPriority
        GetExecutionPriority() const
        {
            return ComponentExecutionPriority::
                Default;
        }

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
