#pragma once

#include <memory>
#include <vector>

#include "PhysiK/API/Handles.h"
#include "PhysiK/Components/CollisionComponent.h"
#include "PhysiK/Math/Vec3.h"

namespace PhysiK
{
    enum class OverlapGeometryType
    {
        Unknown = 0,
        Tetrahedron = 1,
        Triangle = 2,
        Sphere = 3,
        Node = 4
    };

    struct CollisionSphereOverlap
    {
        OverlapGeometryType geometryType = OverlapGeometryType::Unknown;

        ComponentHandle component;
        int primitiveIndex = -1;

        int node0 = -1;
        int node1 = -1;
        int node2 = -1;
        int node3 = -1;

        int overlappedNodeMask = 0;
        int overlappedNodeCount = 0;

        Vec3 sphereCenter;
        float sphereRadius = 0.0f;
        float minDistance = 0.0f;
    };

    class CollisionSphereComponent : public CollisionComponent
    {
    public:
        static std::unique_ptr<CollisionSphereComponent> Create(
            const Vec3& position,
            float radius);

        float radius = 0.5f;

        void SetConnectionSettings(float stiffness, float damping);
        void SetConnectionStiffness(float stiffness);
        float GetConnectionStiffness() const;
        void SetConnectionDamping(float damping);
        float GetConnectionDamping() const;

        void QueryContacts(
            World& world,
            CollisionDetectionEngine& collisionDetectionEngine,
            std::vector<Contact>& outContacts) override;

        void QueryOverlaps(
            World& world,
            std::vector<CollisionSphereOverlap>& outOverlaps) const;
    };
}
