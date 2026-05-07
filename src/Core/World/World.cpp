#include "PhysiK/Core/World/World.h"

#include <algorithm>
#include <cassert>

namespace PhysiK
{
    namespace
    {
        Vec3 WeightedPoint(
            const Node& node0,
            const Node& node1,
            const Node& node2,
            const Node& node3,
            const Vec4& weights)
        {
            return node0.position * weights.x +
                node1.position * weights.y +
                node2.position * weights.z +
                node3.position * weights.w;
        }

        Vec3 WeightedVelocity(
            const Node& node0,
            const Node& node1,
            const Node& node2,
            const Node& node3,
            const Vec4& weights)
        {
            return node0.velocity * weights.x +
                node1.velocity * weights.y +
                node2.velocity * weights.z +
                node3.velocity * weights.w;
        }
    }

    void World::Step(float frameDt)
    {
        if (frameDt <= 0.0f)
        {
            return;
        }

        RunExternalLogic();

        const int steps = std::max(1, substepCount);
        const float substepDt = frameDt / static_cast<float>(steps);

        for (int i = 0; i < steps; ++i)
        {
            ClearForces();
            GenerateCollisionConnections();
            ApplyPointConnectionForces();
            Integrate(substepDt);
            pointConnections.clear();
        }
    }

    int World::AddNode(const Vec3& position, float inverseMass)
    {
        Node node;
        node.position = position;
        node.inverseMass = std::max(0.0f, inverseMass);
        nodes.push_back(node);
        return static_cast<int>(nodes.size()) - 1;
    }

    int World::AddTet(int node0, int node1, int node2, int node3)
    {
        tets.push_back(Tet{node0, node1, node2, node3});
        return static_cast<int>(tets.size()) - 1;
    }

    TetMeshComponent& World::CreateTetMeshComponent(
        const int* nodeIndices,
        int nodeCount,
        const int* tetIndices,
        int tetCount)
    {
        auto component = std::make_unique<TetMeshComponent>();

        if (nodeIndices != nullptr && nodeCount > 0)
        {
            component->nodeIndices.assign(nodeIndices, nodeIndices + nodeCount);
        }

        if (tetIndices != nullptr && tetCount > 0)
        {
            component->tetIndices.assign(tetIndices, tetIndices + tetCount);
        }

        TetMeshComponent& componentRef = *component;
        tetMeshes.push_back(&componentRef);
        components.push_back(std::move(component));
        return componentRef;
    }

    CollisionSphereComponent& World::CreateCollisionSphereComponent(
        const Vec3& position,
        float radius)
    {
        auto component = std::make_unique<CollisionSphereComponent>();
        component->transform.position = position;
        component->radius = std::max(0.0f, radius);

        CollisionSphereComponent& componentRef = *component;
        collisionComponents.push_back(&componentRef);
        components.push_back(std::move(component));
        return componentRef;
    }

    void World::AddPointConnection(const PointConnection& connection)
    {
        if (HasValidNodeIndices(connection))
        {
            pointConnections.push_back(connection);
        }
    }

    void World::SetExternalLogicCallback(ExternalLogicCallback callback, void* userData)
    {
        externalLogicCallback = callback;
        externalLogicUserData = userData;
    }

    void World::ClearExternalLogicCallback()
    {
        externalLogicCallback = nullptr;
        externalLogicUserData = nullptr;
    }

    void World::SetSubstepCount(int count)
    {
        substepCount = std::max(1, count);
    }

    int World::GetSubstepCount() const
    {
        return substepCount;
    }

    Node& World::GetNode(int index)
    {
        assert(index >= 0 && index < static_cast<int>(nodes.size()));
        return nodes[static_cast<std::size_t>(index)];
    }

    const Node& World::GetNode(int index) const
    {
        assert(index >= 0 && index < static_cast<int>(nodes.size()));
        return nodes[static_cast<std::size_t>(index)];
    }

    const std::vector<Tet>& World::GetTets() const
    {
        return tets;
    }

    const std::vector<PointConnection>& World::GetPointConnections() const
    {
        return pointConnections;
    }

    void World::RunExternalLogic()
    {
        if (externalLogicCallback != nullptr)
        {
            externalLogicCallback(static_cast<WorldHandle>(this), externalLogicUserData);
        }
    }

    void World::GenerateCollisionConnections()
    {
        std::vector<Contact> contacts;

        for (CollisionComponent* component : collisionComponents)
        {
            if (component == nullptr || component->isSensor || !component->generateConnections)
            {
                continue;
            }

            contacts.clear();
            collisionDetectionEngine.QueryContacts(*this, *component, contacts);

            for (const Contact& contact : contacts)
            {
                AddPointConnectionFromContact(contact);
            }
        }
    }

    void World::AddPointConnectionFromContact(const Contact& contact)
    {
        if (contact.penetrationDepth <= 0.0f)
        {
            return;
        }

        PointConnection connection;
        connection.node0 = contact.tetNode0;
        connection.node1 = contact.tetNode1;
        connection.node2 = contact.tetNode2;
        connection.node3 = contact.tetNode3;
        connection.barycentric = contact.barycentric;
        connection.targetPosition = contact.worldPoint + contact.normal * contact.penetrationDepth;
        connection.stiffness = contact.stiffness;
        connection.damping = contact.damping;
        AddPointConnection(connection);
    }

    void World::ApplyPointConnectionForces()
    {
        for (const PointConnection& connection : pointConnections)
        {
            Node& node0 = nodes[static_cast<std::size_t>(connection.node0)];
            Node& node1 = nodes[static_cast<std::size_t>(connection.node1)];
            Node& node2 = nodes[static_cast<std::size_t>(connection.node2)];
            Node& node3 = nodes[static_cast<std::size_t>(connection.node3)];

            const Vec3 point = WeightedPoint(node0, node1, node2, node3, connection.barycentric);
            const Vec3 velocity = WeightedVelocity(node0, node1, node2, node3, connection.barycentric);
            const Vec3 springForce = (connection.targetPosition - point) * connection.stiffness;
            const Vec3 dampingForce = velocity * (-connection.damping);
            const Vec3 pointForce = springForce + dampingForce;

            node0.force += pointForce * connection.barycentric.x;
            node1.force += pointForce * connection.barycentric.y;
            node2.force += pointForce * connection.barycentric.z;
            node3.force += pointForce * connection.barycentric.w;
        }
    }

    void World::Integrate(float dt)
    {
        for (Node& node : nodes)
        {
            if (node.inverseMass <= 0.0f)
            {
                continue;
            }

            const Vec3 acceleration = node.force * node.inverseMass;
            node.velocity += acceleration * dt;
            node.position += node.velocity * dt;
        }
    }

    void World::ClearForces()
    {
        for (Node& node : nodes)
        {
            node.force = Vec3{};
        }
    }

    bool World::HasValidNodeIndices(const PointConnection& connection) const
    {
        const int nodeCount = static_cast<int>(nodes.size());
        return connection.node0 >= 0 && connection.node0 < nodeCount &&
            connection.node1 >= 0 && connection.node1 < nodeCount &&
            connection.node2 >= 0 && connection.node2 < nodeCount &&
            connection.node3 >= 0 && connection.node3 < nodeCount;
    }
}
