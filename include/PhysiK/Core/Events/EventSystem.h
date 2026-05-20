#pragma once

#include <vector>

#include "PhysiK/API/PhysiKAPI.h"
#include "PhysiK/Core/Events/PhysicsEvent.h"

namespace PhysiK
{
    class Component;

    class PHYSIK_API EventSystem
    {
    public:
        void Subscribe(Component* listener, PhysicsEventType type);
        void Unsubscribe(Component* listener, PhysicsEventType type);
        void UnsubscribeAll(Component* listener);
        void Emit(const PhysicsEvent& event);

    private:
        struct Subscription
        {
            Component* listener = nullptr;
            PhysicsEventType type = PhysicsEventType::TetMeshTopologyChanged;
        };

        std::vector<Subscription> subscriptions;
    };
}
