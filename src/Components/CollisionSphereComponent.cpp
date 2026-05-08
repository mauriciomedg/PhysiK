#include "PhysiK/Components/CollisionSphereComponent.h"

#include <cmath>

#include "PhysiK/Components/TetMeshComponent.h"
#include "PhysiK/Core/Collision/CollisionDetectionEngine.h"
#include "PhysiK/Core/World/World.h"

namespace PhysiK
{
    namespace
    {
        Vec3 NormalizeOrFallback(const Vec3& value)
        {
            const float length = value.Length();
            if (length > 0.00001f)
            {
                return value / length;
            }

            return Vec3{0.0f, 0.0f, 1.0f};
        }
    }

    void CollisionSphereComponent::QueryContacts(
        World& world,
        CollisionDetectionEngine& collisionDetectionEngine,
        std::vector<Contact>& outContacts)
    {
        (void)collisionDetectionEngine;

        if (!active || radius <= 0.0f)
        {
            return;
        }

        const Vec4 centroidWeights{0.25f, 0.25f, 0.25f, 0.25f};

        for (const TetMeshComponent* tetMesh : world.GetTetMeshes())
        {
            if (tetMesh == nullptr || !tetMesh->active)
            {
                continue;
            }

            for (const Tet& tet : tetMesh->tets)
            {
                const Node& node0 = world.GetNode(tet.node0);
                const Node& node1 = world.GetNode(tet.node1);
                const Node& node2 = world.GetNode(tet.node2);
                const Node& node3 = world.GetNode(tet.node3);

                const Vec3 point = node0.position * centroidWeights.x +
                    node1.position * centroidWeights.y +
                    node2.position * centroidWeights.z +
                    node3.position * centroidWeights.w;
                const Vec3 centerToPoint = point - transform.position;
                const float distance = centerToPoint.Length();

                if (distance >= radius)
                {
                    continue;
                }

                const Vec3 normal = NormalizeOrFallback(centerToPoint);

                Contact contact;
                contact.node0 = tet.node0;
                contact.node1 = tet.node1;
                contact.node2 = tet.node2;
                contact.node3 = tet.node3;
                contact.barycentric = centroidWeights;
                contact.worldPoint = point;
                contact.normal = normal;
                contact.penetrationDepth = radius - distance;
                contact.stiffness = contactStiffness;
                contact.damping = contactDamping;
                outContacts.push_back(contact);
            }
        }
    }
}
