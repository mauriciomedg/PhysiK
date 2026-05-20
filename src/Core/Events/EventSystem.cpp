#include "PhysiK/Core/Events/EventSystem.h"

#include <algorithm>

#include "PhysiK/Components/Component.h"

namespace PhysiK
{
    void EventSystem::Subscribe(Component* listener, PhysicsEventType type)
    {
        if (listener == nullptr)
        {
            return;
        }

        const auto existing = std::find_if(
            subscriptions.begin(),
            subscriptions.end(),
            [listener, type](const Subscription& subscription)
            {
                return subscription.listener == listener && subscription.type == type;
            });
        if (existing != subscriptions.end())
        {
            return;
        }

        subscriptions.push_back(Subscription{listener, type});
    }

    void EventSystem::Unsubscribe(Component* listener, PhysicsEventType type)
    {
        subscriptions.erase(
            std::remove_if(
                subscriptions.begin(),
                subscriptions.end(),
                [listener, type](const Subscription& subscription)
                {
                    return subscription.listener == listener && subscription.type == type;
                }),
            subscriptions.end());
    }

    void EventSystem::UnsubscribeAll(Component* listener)
    {
        subscriptions.erase(
            std::remove_if(
                subscriptions.begin(),
                subscriptions.end(),
                [listener](const Subscription& subscription)
                {
                    return subscription.listener == listener;
                }),
            subscriptions.end());
    }

    void EventSystem::Emit(const PhysicsEvent& event)
    {
        for (const Subscription& subscription : subscriptions)
        {
            if (subscription.listener != nullptr && subscription.type == event.type)
            {
                subscription.listener->OnPhysicsEvent(event);
            }
        }
    }
}
