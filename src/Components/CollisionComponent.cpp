#include "PhysiK/Components/CollisionComponent.h"

namespace PhysiK
{
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
