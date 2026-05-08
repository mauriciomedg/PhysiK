#include "PhysiK/Components/CollisionComponent.h"

namespace PhysiK
{
    void CollisionComponent::QueryContacts(
        World& world,
        CollisionDetectionEngine& collisionDetectionEngine,
        std::vector<Contact>& outContacts)
    {
        (void)world;
        (void)collisionDetectionEngine;
        (void)outContacts;
    }

    void CollisionComponent::UpdateKinematicTarget(World& world)
    {
        (void)world;

        Transform target;
        if (ConsumeKinematicTarget(target))
        {
            transform = target;
        }
    }

    void CollisionComponent::SetKinematicTarget(const Transform& target)
    {
        kinematicTarget = target;
        hasKinematicTarget = true;
    }

    bool CollisionComponent::ConsumeKinematicTarget(Transform& outTarget)
    {
        if (!hasKinematicTarget)
        {
            return false;
        }

        outTarget = kinematicTarget;
        hasKinematicTarget = false;
        return true;
    }
}
