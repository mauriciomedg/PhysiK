#pragma once

#include "PhysiK/Components/CollisionComponent.h"

namespace PhysiK
{
    class CollisionSphereComponent : public CollisionComponent
    {
    public:
        float radius = 0.5f;

        void QueryContacts(
            World& world,
            CollisionDetectionEngine& collisionDetectionEngine,
            std::vector<Contact>& outContacts) override;
    };
}
