#include "PhysiK/Components/CollisionSphereComponent.h"

#include <algorithm>
#include <cmath>
#include <memory>

#include "PhysiK/Components/TetMeshComponent.h"
#include "PhysiK/Components/TetMeshPhysicsComponent.h"
#include "PhysiK/Core/PhysicsConnections/PointConnection.h"
#include "PhysiK/Core/Solvers/SolverData.h"
#include "PhysiK/Core/World/World.h"
#include "PhysiK/PhysicsData/Contact.h"

namespace PhysiK
{
    namespace
    {
        bool HasValidWorldNodes(const int nodes[4], const World& world)
        {
            const int nodeCount = static_cast<int>(world.GetNodes().size());
            return nodes[0] >= 0 && nodes[0] < nodeCount &&
                nodes[1] >= 0 && nodes[1] < nodeCount &&
                nodes[2] >= 0 && nodes[2] < nodeCount &&
                nodes[3] >= 0 && nodes[3] < nodeCount;
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
            const auto* tetMesh =
                dynamic_cast<const TetMeshPhysicsComponent*>(component.get());
            if (tetMesh == nullptr || !tetMesh->active)
            {
                continue;
            }

            for (int tetIndex = 0; tetIndex < tetMesh->GetTetCount(); ++tetIndex)
            {
                if (!tetMesh->IsTetActive(tetIndex))
                {
                    continue;
                }

                const int localNodes[4] = {
                    tetMesh->GetTetNodeIndex(tetIndex, 0),
                    tetMesh->GetTetNodeIndex(tetIndex, 1),
                    tetMesh->GetTetNodeIndex(tetIndex, 2),
                    tetMesh->GetTetNodeIndex(tetIndex, 3)};
                const int worldNodes[4] = {
                    tetMesh->GetWorldNodeIndex(localNodes[0]),
                    tetMesh->GetWorldNodeIndex(localNodes[1]),
                    tetMesh->GetWorldNodeIndex(localNodes[2]),
                    tetMesh->GetWorldNodeIndex(localNodes[3])};

                if (!HasValidWorldNodes(worldNodes, world))
                {
                    continue;
                }

                const Node& node0 = world.GetNode(worldNodes[0]);
                const Node& node1 = world.GetNode(worldNodes[1]);
                const Node& node2 = world.GetNode(worldNodes[2]);
                const Node& node3 = world.GetNode(worldNodes[3]);

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
                contact.node0 = worldNodes[0];
                contact.node1 = worldNodes[1];
                contact.node2 = worldNodes[2];
                contact.node3 = worldNodes[3];
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
            const auto* tetMesh =
                dynamic_cast<const TetMeshPhysicsComponent*>(component.get());
            if (tetMesh == nullptr || !tetMesh->active)
            {
                continue;
            }

            for (int tetIndex = 0; tetIndex < tetMesh->GetTetCount(); ++tetIndex)
            {
                if (!tetMesh->IsTetActive(tetIndex))
                {
                    continue;
                }

                const int localNodes[4] = {
                    tetMesh->GetTetNodeIndex(tetIndex, 0),
                    tetMesh->GetTetNodeIndex(tetIndex, 1),
                    tetMesh->GetTetNodeIndex(tetIndex, 2),
                    tetMesh->GetTetNodeIndex(tetIndex, 3)};
                const int nodes[4] = {
                    tetMesh->GetWorldNodeIndex(localNodes[0]),
                    tetMesh->GetWorldNodeIndex(localNodes[1]),
                    tetMesh->GetWorldNodeIndex(localNodes[2]),
                    tetMesh->GetWorldNodeIndex(localNodes[3])};
                if (!HasValidWorldNodes(nodes, world))
                {
                    continue;
                }
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
                overlap.node0 = nodes[0];
                overlap.node1 = nodes[1];
                overlap.node2 = nodes[2];
                overlap.node3 = nodes[3];
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
