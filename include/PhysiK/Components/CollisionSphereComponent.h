#pragma once

#include <memory>

#include "PhysiK/Components/CollisionComponent.h"
#include "PhysiK/Math/Vec3.h"

namespace PhysiK
{
    class CollisionSphereComponent : public CollisionComponent
    {
    public:
        static std::unique_ptr<CollisionSphereComponent> Create(
            const Vec3& position,
            float radius);

        float radius = 0.5f;

        void QueryContacts(
            World& world,
            CollisionDetectionEngine& collisionDetectionEngine,
            std::vector<Contact>& outContacts) override;
    };
}
