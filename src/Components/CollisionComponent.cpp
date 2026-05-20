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

    void CollisionComponent::PreUpdate(World& world, float dt)
    {
        (void)world;
        (void)dt;

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
