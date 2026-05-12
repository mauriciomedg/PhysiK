#include "PhysiK/Core/World/World.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <memory>

namespace PhysiK
{
    namespace
    {
        bool IsFinite(float value)
        {
            return std::isfinite(value);
        }

        bool IsFinite(const Vec3& value)
        {
            return IsFinite(value.x) && IsFinite(value.y) && IsFinite(value.z);
        }

    }

    World::World() = default;

    void World::Step(float frameDt)
    {
        if (frameDt <= 0.0f)
        {
            return;
        }

        RunExternalLogic();
        UpdateFrameComponents(frameDt);
        UpdateKinematicTargets();

        const int steps = std::max(1, substepCount);
        const float substepDt = frameDt / static_cast<float>(steps);

        for (int i = 0; i < steps; ++i)
        {
            SolverData solverData;
            PrecomputeSolve(solverData, substepDt);

            if (solverMode == SolverMode::ImplicitEuler)
            {
                if (SolveImplicitLinearSystem(solverData, substepDt))
                {
                    IntegrateImplicitEuler(solverData, substepDt);
                }
            }
            else
            {
                IntegrateExplicitEuler(solverData, substepDt);
            }

            ClearTransientConnections();
        }
    }

    int World::AddNode(const Vec3& position)
    {
        Node node;
        node.position = position;
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

    void World::SetSolverMode(SolverMode mode)
    {
        solverMode = mode;
    }

    SolverMode World::GetSolverMode() const
    {
        return solverMode;
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

    void World::SetNodeFixed(int nodeIndex, bool fixed)
    {
        assert(nodeIndex >= 0 && nodeIndex < static_cast<int>(nodes.size()));
        Node& node = nodes[static_cast<std::size_t>(nodeIndex)];
        if (fixed)
        {
            node.fixed = true;
            node.velocity = Vec3{};
            return;
        }

        node.fixed = false;
    }

    bool World::IsNodeFixed(int nodeIndex) const
    {
        assert(nodeIndex >= 0 && nodeIndex < static_cast<int>(nodes.size()));
        return nodes[static_cast<std::size_t>(nodeIndex)].fixed;
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

    void World::RunExternalLogic()
    {
        if (externalLogicCallback != nullptr)
        {
            externalLogicCallback(static_cast<WorldHandle>(this), externalLogicUserData);
        }
    }

    void World::UpdateFrameComponents(float frameDt)
    {
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
            if (component != nullptr && component->active)
            {
                component->UpdateKinematicTarget(*this);
            }
        }
    }

    void World::BuildSolverData(SolverData& solverData, float dt)
    {
        GenerateCollisionConnections();

        AssembleComponentSystems(solverData, dt);
        AddDefaultNodeMasses(solverData);
        AddGravityForces(solverData);
        AssembleConnectionSystems(solverData, dt);
    }

    void World::AddDefaultNodeMasses(SolverData& solverData)
    {
        solverData.AssembleMasses(static_cast<int>(nodes.size()));

        for (int i = 0; i < static_cast<int>(nodes.size()); ++i)
        {
            if (nodes[static_cast<std::size_t>(i)].fixed ||
                solverData.HasNodeMassContribution(i) ||
                solverData.GetAssembledMassForNode(i) > 0.0f)
            {
                continue;
            }

            solverData.AddNodeMass(i, 1.0f);
        }
    }

    void World::AddGravityForces(SolverData& solverData)
    {
        solverData.AssembleMasses(static_cast<int>(nodes.size()));

        for (int i = 0; i < static_cast<int>(nodes.size()); ++i)
        {
            if (nodes[static_cast<std::size_t>(i)].fixed)
            {
                continue;
            }

            const float mass = solverData.GetAssembledMassForNode(i);
            if (!std::isfinite(mass) || mass <= 0.0f)
            {
                continue;
            }

            solverData.AddNodeForce(i, gravity * mass);
        }
    }

    void World::AssembleConnectionSystems(SolverData& solverData, float dt)
    {
        for (const std::unique_ptr<PhysicsConnection>& connection : transientConnections)
        {
            if (connection != nullptr)
            {
                connection->UpdateSystem(*this, solverData, dt);
            }
        }
    }

    void World::GenerateCollisionConnections()
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
                GeneratePointConnectionFromContact(contact);
            }
        }
    }

    void World::AssembleComponentSystems(SolverData& solverData, float dt)
    {
        for (const std::unique_ptr<Component>& component : components)
        {
            if (component != nullptr && component->active)
            {
                component->UpdateSystem(*this, solverData, dt);
            }
        }
    }

    void World::GeneratePointConnectionFromContact(const Contact& contact)
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
    }

    void World::PrecomputeSolve(SolverData& solverData, float dt)
    {
        solverData.Clear();
        BuildSolverData(solverData, dt);
        solverData.PrecomputeImplicitSolve(nodes, dt);
    }

    bool World::SolveImplicitLinearSystem(SolverData& solverData, float dt)
    {
        (void)dt;
        return solverData.SolveImplicitLinearSystem();
    }

    bool World::IntegrateImplicitEuler(const SolverData& solverData, float dt)
    {
        std::vector<Vec3> updatedVelocities(nodes.size());
        std::vector<Vec3> updatedPositions(nodes.size());

        for (int nodeIndex = 0; nodeIndex < static_cast<int>(nodes.size()); ++nodeIndex)
        {
            const int dynamicBlock = solverData.GetDynamicBlockForNode(nodeIndex);
            if (dynamicBlock < 0)
            {
                continue;
            }

            const int baseDof = dynamicBlock * 3;
            const Node& node = nodes[static_cast<std::size_t>(nodeIndex)];
            const std::vector<float>& deltaVelocity = solverData.GetDeltaVelocity();
            if (baseDof + 2 >= static_cast<int>(deltaVelocity.size()))
            {
                return false;
            }

            const Vec3 velocityChange{
                deltaVelocity[static_cast<std::size_t>(baseDof + 0)],
                deltaVelocity[static_cast<std::size_t>(baseDof + 1)],
                deltaVelocity[static_cast<std::size_t>(baseDof + 2)]};
            if (!IsFinite(node.position) || !IsFinite(node.velocity) || !IsFinite(velocityChange))
            {
                return false;
            }

            const Vec3 updatedVelocity = node.velocity + velocityChange;
            const Vec3 updatedPosition = node.position + updatedVelocity * dt;
            if (!IsFinite(updatedVelocity) || !IsFinite(updatedPosition))
            {
                return false;
            }

            updatedVelocities[static_cast<std::size_t>(nodeIndex)] = updatedVelocity;
            updatedPositions[static_cast<std::size_t>(nodeIndex)] = updatedPosition;
        }

        for (int nodeIndex = 0; nodeIndex < static_cast<int>(nodes.size()); ++nodeIndex)
        {
            if (solverData.GetDynamicBlockForNode(nodeIndex) < 0)
            {
                continue;
            }

            Node& node = nodes[static_cast<std::size_t>(nodeIndex)];
            node.velocity = updatedVelocities[static_cast<std::size_t>(nodeIndex)];
            node.position = updatedPositions[static_cast<std::size_t>(nodeIndex)];
        }

        return true;
    }

    void World::IntegrateExplicitEuler(const SolverData& solverData, float dt)
    {
        const std::vector<float>& masses = solverData.GetAssembledMasses();
        const std::vector<Vec3>& forces = solverData.GetAssembledForces();

        for (int i = 0; i < static_cast<int>(nodes.size()); ++i)
        {
            Node& node = nodes[static_cast<std::size_t>(i)];
            if (node.fixed ||
                i >= static_cast<int>(masses.size()) ||
                i >= static_cast<int>(forces.size()))
            {
                continue;
            }

            const float mass = masses[static_cast<std::size_t>(i)];
            if (!std::isfinite(mass) || mass <= 0.0f)
            {
                continue;
            }

            const Vec3 acceleration = forces[static_cast<std::size_t>(i)] / mass;
            node.velocity += acceleration * dt;
            node.position += node.velocity * dt;
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
