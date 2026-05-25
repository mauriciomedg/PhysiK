#include "PhysiK/Components/CollisionSphereComponent.h"

#include <algorithm>
#include <cmath>
#include <memory>

#include "PhysiK/Components/TetMeshComponent.h"
#include "PhysiK/Core/PhysicsConnections/PointConnection.h"
#include "PhysiK/Core/Solvers/SolverData.h"
#include "PhysiK/Core/World/World.h"
#include "PhysiK/PhysicsData/Contact.h"

namespace PhysiK
{
    namespace
    {
        bool HasValidTetNodes(const Tet& tet, const World& world)
        {
            const int nodeCount = static_cast<int>(world.GetNodes().size());
            return tet.node0 >= 0 && tet.node0 < nodeCount &&
                tet.node1 >= 0 && tet.node1 < nodeCount &&
                tet.node2 >= 0 && tet.node2 < nodeCount &&
                tet.node3 >= 0 && tet.node3 < nodeCount;
        }

        Vec3 NormalizeOrFallback(const Vec3& value)
        {
            const float length = value.Length();
            if (length > 0.00001f)
            {
                return value / length;
            }

            return Vec3{0.0f, 0.0f, 1.0f};
        }

        float SanitizeConnectionValue(float value)
        {
            if (!std::isfinite(value))
            {
                return 0.0f;
            }

            return std::max(0.0f, value);
        }

        void AddPointConnectionFromContact(World& world, const Contact& contact)
        {
            if (contact.penetrationDepth <= 0.0f)
            {
                return;
            }

            PointConnection connection;
            connection.node0 = contact.node0;
            connection.node1 = contact.node1;
            connection.node2 = contact.node2;
            connection.node3 = contact.node3;
            connection.barycentric = contact.barycentric;
            connection.targetPosition =
                contact.worldPoint + contact.normal * contact.penetrationDepth;
            connection.stiffness = contact.stiffness;
            connection.damping = contact.damping;
            world.AddPointConnection(connection);
        }
    }

    std::unique_ptr<CollisionSphereComponent> CollisionSphereComponent::Create(
        const Vec3& position,
        float radius)
    {
        auto component = std::make_unique<CollisionSphereComponent>();
        component->transform.position = position;
        component->radius = std::max(0.0f, radius);
        return component;
    }

    void CollisionSphereComponent::SetConnectionSettings(float stiffness, float damping)
    {
        SetConnectionStiffness(stiffness);
        SetConnectionDamping(damping);
    }

    void CollisionSphereComponent::SetConnectionStiffness(float stiffness)
    {
        contactStiffness = SanitizeConnectionValue(stiffness);
    }

    float CollisionSphereComponent::GetConnectionStiffness() const
    {
        return contactStiffness;
    }

    void CollisionSphereComponent::SetConnectionDamping(float damping)
    {
        contactDamping = SanitizeConnectionValue(damping);
    }

    float CollisionSphereComponent::GetConnectionDamping() const
    {
        return contactDamping;
    }

    void CollisionSphereComponent::UpdateSystem(
        World& world,
        SolverData& solverData,
        float dt)
    {
        (void)solverData;
        (void)dt;

        if (!active || radius <= 0.0f)
        {
            return;
        }

        const Vec4 centroidWeights{0.25f, 0.25f, 0.25f, 0.25f};

        if (isSensor || !generateConnections)
        {
            return;
        }

        for (const std::unique_ptr<Component>& component : world.GetComponents())
        {
            const auto* tetMesh = dynamic_cast<const TetMeshComponent*>(component.get());
            if (tetMesh == nullptr || !tetMesh->active)
            {
                continue;
            }

            for (const Tet& tet : tetMesh->tets)
            {
                if (!tet.active)
                {
                    continue;
                }

                if (!HasValidTetNodes(tet, world))
                {
                    continue;
                }

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
                contact.stiffness = GetConnectionStiffness();
                contact.damping = GetConnectionDamping();
                AddPointConnectionFromContact(world, contact);
            }
        }
    }

    void CollisionSphereComponent::QueryOverlaps(
        World& world,
        std::vector<CollisionSphereOverlap>& outOverlaps) const
    {
        if (!active || radius <= 0.0f)
        {
            return;
        }

        const Vec3 sphereCenter = transform.position;
        const float radiusSquared = radius * radius;
        const std::vector<std::unique_ptr<Component>>& components = world.GetComponents();

        for (int componentIndex = 0; componentIndex < static_cast<int>(components.size()); ++componentIndex)
        {
            const std::unique_ptr<Component>& component =
                components[static_cast<std::size_t>(componentIndex)];
            const auto* tetMesh = dynamic_cast<const TetMeshComponent*>(component.get());
            if (tetMesh == nullptr || !tetMesh->active)
            {
                continue;
            }

            for (int tetIndex = 0; tetIndex < static_cast<int>(tetMesh->tets.size()); ++tetIndex)
            {
                const Tet& tet = tetMesh->tets[static_cast<std::size_t>(tetIndex)];
                if (!tet.active || !HasValidTetNodes(tet, world))
                {
                    continue;
                }

                const int nodes[4] = {tet.node0, tet.node1, tet.node2, tet.node3};
                int overlappedNodeMask = 0;
                int overlappedNodeCount = 0;
                float minNodeDistance = 0.0f;
                bool hasDistance = false;

                for (int localNode = 0; localNode < 4; ++localNode)
                {
                    const Vec3 difference =
                        world.GetNode(nodes[localNode]).position - sphereCenter;
                    const float distanceSquared = difference.LengthSquared();
                    const float distance = std::sqrt(std::max(0.0f, distanceSquared));
                    if (!hasDistance || distance < minNodeDistance)
                    {
                        minNodeDistance = distance;
                        hasDistance = true;
                    }

                    if (distanceSquared <= radiusSquared)
                    {
                        overlappedNodeMask |= (1 << localNode);
                        ++overlappedNodeCount;
                    }
                }

                if (overlappedNodeCount <= 0)
                {
                    continue;
                }

                CollisionSphereOverlap overlap;
                overlap.geometryType = OverlapGeometryType::Tetrahedron;
                overlap.component = world.GetComponentHandleByIndex(componentIndex);
                overlap.primitiveIndex = tetIndex;
                overlap.node0 = tet.node0;
                overlap.node1 = tet.node1;
                overlap.node2 = tet.node2;
                overlap.node3 = tet.node3;
                overlap.overlappedNodeMask = overlappedNodeMask;
                overlap.overlappedNodeCount = overlappedNodeCount;
                overlap.sphereCenter = sphereCenter;
                overlap.sphereRadius = radius;
                overlap.minDistance = minNodeDistance;
                outOverlaps.push_back(overlap);
            }
        }
    }
}
