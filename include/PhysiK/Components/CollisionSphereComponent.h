#pragma once

#include <memory>
#include <vector>

#include "PhysiK/API/Handles.h"
#include "PhysiK/Components/CollisionComponent.h"
#include "PhysiK/Math/Vec3.h"

namespace PhysiK
{
    struct SphereTetOverlap
    {
        ComponentHandle tetMeshComponent;
        int tetIndex = -1;

        int node0 = -1;
        int node1 = -1;
        int node2 = -1;
        int node3 = -1;

        int overlappedNodeMask = 0;
        int overlappedNodeCount = 0;

        Vec3 sphereCenter;
        float sphereRadius = 0.0f;
        float minNodeDistance = 0.0f;
    };

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

        void QueryOverlappingTets(
            World& world,
            std::vector<SphereTetOverlap>& outOverlaps) const;
    };
}
