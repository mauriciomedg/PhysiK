#include "PhysiK/Core/World/World.h"

#include <algorithm>
#include <cassert>
#include <memory>

#include "PhysiK/Components/CollisionComponent.h"

namespace PhysiK
{
    World::World() = default;

    void World::Step(float frameDt)
    {
        if (frameDt <= 0.0f)
        {
            return;
        }

        RunExternalLogic(frameDt);
        UpdateKinematicTargets();

        const int steps = std::max(1, substepCount);
        const float substepDt = frameDt / static_cast<float>(steps);

        for (int i = 0; i < steps; ++i)
        {
            AccumulateForces(substepDt);
            Integrate(substepDt);
            ClearTransientConnections();
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

    Component* World::GetComponent(ComponentHandle handle)
    {
        if (!IsComponentHandleValid(handle))
        {
            return nullptr;
        }

        return components[handle.index].get();
    }

    const Component* World::GetComponent(ComponentHandle handle) const
    {
        if (!IsComponentHandleValid(handle))
        {
            return nullptr;
        }

        return components[handle.index].get();
    }

    void World::DestroyComponent(ComponentHandle handle)
    {
        if (!IsComponentHandleValid(handle))
        {
            return;
        }

        components[handle.index].reset();
        ++componentGenerations[handle.index];

        if (componentGenerations[handle.index] == 0u)
        {
            componentGenerations[handle.index] = 1u;
        }

        freeComponentSlots.push_back(handle.index);
    }

    void World::AddPointConnection(const PointConnection& connection)
    {
        if (HasValidNodeIndices(connection))
        {
            transientConnections.push_back(std::make_unique<PointConnection>(connection));
        }
    }

    void World::AddTransientConnection(std::unique_ptr<PhysicsConnection> connection)
    {
        if (connection != nullptr)
        {
            transientConnections.push_back(std::move(connection));
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

    void World::SetGravity(const Vec3& value)
    {
        gravity = value;
    }

    const Vec3& World::GetGravity() const
    {
        return gravity;
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

    const std::vector<Node>& World::GetNodes() const
    {
        return nodes;
    }

    void World::SetNodePosition(int index, const Vec3& position)
    {
        assert(index >= 0 && index < static_cast<int>(nodes.size()));
        Node& node = nodes[static_cast<std::size_t>(index)];
        node.position = position;
        node.velocity = Vec3{};
    }

    const std::vector<std::unique_ptr<Component>>& World::GetComponents() const
    {
        return components;
    }

    int World::GetTransientConnectionCount() const
    {
        return static_cast<int>(transientConnections.size());
    }

    ComponentHandle World::AddComponent(std::unique_ptr<Component> component)
    {
        if (component == nullptr)
        {
            return ComponentHandle{};
        }

        if (!freeComponentSlots.empty())
        {
            const std::uint32_t slotIndex = freeComponentSlots.back();
            freeComponentSlots.pop_back();
            components[slotIndex] = std::move(component);
            return ComponentHandle{slotIndex, componentGenerations[slotIndex]};
        }

        components.push_back(std::move(component));
        componentGenerations.push_back(1u);

        return ComponentHandle{
            static_cast<std::uint32_t>(components.size() - 1u),
            componentGenerations.back()};
    }

    bool World::IsComponentHandleValid(ComponentHandle handle) const
    {
        if (!handle.IsValid() || handle.index >= components.size())
        {
            return false;
        }

        return componentGenerations[handle.index] == handle.generation &&
            components[handle.index] != nullptr;
    }

    void World::RunExternalLogic(float frameDt)
    {
        if (externalLogicCallback != nullptr)
        {
            externalLogicCallback(static_cast<WorldHandle>(this), externalLogicUserData);
        }

        for (const std::unique_ptr<Component>& component : components)
        {
            if (component != nullptr && component->active)
            {
                component->UpdateFrame(*this, frameDt);
            }
        }
    }

    void World::UpdateKinematicTargets()
    {
        for (const std::unique_ptr<Component>& component : components)
        {
            if (component == nullptr)
            {
                continue;
            }

            if (auto* collision = dynamic_cast<CollisionComponent*>(component.get()))
            {
                Transform target;
                if (collision->ConsumeKinematicTarget(target))
                {
                    collision->transform = target;
                }
            }
        }
    }

    void World::AccumulateForces(float dt)
    {
        ClearForces();
        SolverData solverData;
        solverData.Clear();
        AddGravityForces(solverData);
        AddCollisionForces(solverData, dt);
        AddPhysicsModelForces(solverData, dt);
        AddConnectionForces(solverData, dt);
        Solve(solverData, dt);
    }

    void World::AddGravityForces(SolverData& solverData)
    {
        for (int i = 0; i < static_cast<int>(nodes.size()); ++i)
        {
            const Node& node = nodes[static_cast<std::size_t>(i)];
            if (node.inverseMass <= 0.0f)
            {
                continue;
            }

            solverData.AddNodeForce(i, gravity / node.inverseMass);
        }
    }

    void World::AddConnectionForces(SolverData& solverData, float dt)
    {
        for (const std::unique_ptr<PhysicsConnection>& connection : transientConnections)
        {
            if (connection != nullptr)
            {
                connection->UpdateSystem(*this, solverData, dt);
            }
        }
    }

    void World::AddCollisionForces(SolverData& solverData, float dt)
    {
        std::vector<Contact> contacts;

        for (const std::unique_ptr<Component>& component : components)
        {
            if (component == nullptr || !component->active)
            {
                continue;
            }

            contacts.clear();
            component->QueryContacts(*this, collisionDetectionEngine, contacts);

            for (const Contact& contact : contacts)
            {
                AddPointConnectionFromContact(contact, solverData, dt);
            }
        }
    }

    void World::AddPhysicsModelForces(SolverData& solverData, float dt)
    {
        for (const std::unique_ptr<Component>& component : components)
        {
            if (component != nullptr && component->active)
            {
                component->UpdateSystem(*this, solverData, dt);
            }
        }
    }

    void World::AddPointConnectionFromContact(
        const Contact& contact,
        SolverData& solverData,
        float dt)
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
        connection.targetPosition = contact.worldPoint + contact.normal * contact.penetrationDepth;
        connection.stiffness = contact.stiffness;
        connection.damping = contact.damping;
        AddPointConnection(connection);
        (void)solverData;
        (void)dt;
    }

    void World::Solve(SolverData& solverData, float dt)
    {
        (void)dt;

        for (const SolverData::NodeForce& nodeForce : solverData.GetNodeForces())
        {
            if (nodeForce.node >= 0 && nodeForce.node < static_cast<int>(nodes.size()))
            {
                nodes[static_cast<std::size_t>(nodeForce.node)].force += nodeForce.force;
            }
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

    void World::ClearTransientConnections()
    {
        transientConnections.clear();
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
