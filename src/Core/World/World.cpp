#include "PhysiK/Core/World/World.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>

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

        constexpr float MinConjugateGradientTolerance = 1.0e-8f;
        constexpr float MaxConjugateGradientTolerance = 1.0e-1f;
        constexpr int MinConjugateGradientMaxIterations = 1;
        constexpr int MaxConjugateGradientMaxIterations = 1024;
    }

    World::World() = default;

    void World::Step(float frameDt)
    {
        if (frameDt < 0.0f)
        {
            return;
        }

        PreUpdateComponents(frameDt);
        MarkSubstepConnectionBegin();
        if (frameDt == 0.0f)
        {
            PostUpdateComponents(frameDt);
            ClearFrameConnections();
            return;
        }

        const int steps = std::max(1, substepCount);
        const float substepDt = frameDt / static_cast<float>(steps);

        SolverData solverData;
        for (int i = 0; i < steps; ++i)
        {
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

            ClearSubstepConnections();
        }

        PostUpdateComponents(frameDt);
        ClearFrameConnections();
    }

    int World::AddNode(const Vec3& position)
    {
        const int stateIndex = state.AddNodeState(position);
        Node node;
        node.stateIndex = stateIndex;
        nodes.push_back(node);
        assert(stateIndex == static_cast<int>(nodes.size()) - 1);
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

        Component* component = components[handle.index].get();
        UnregisterComponentFromExecution(component);
        eventSystem.UnsubscribeAll(component);
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

    void World::SubscribeToEvent(Component* listener, PhysicsEventType type)
    {
        eventSystem.Subscribe(listener, type);
    }

    void World::UnsubscribeFromEvent(Component* listener, PhysicsEventType type)
    {
        eventSystem.Unsubscribe(listener, type);
    }

    void World::EmitEvent(const PhysicsEvent& event)
    {
        eventSystem.Emit(event);
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

    void World::SetConjugateGradientTolerance(float tolerance)
    {
        if (!std::isfinite(tolerance))
        {
            return;
        }

        conjugateGradientSettings.tolerance = std::clamp(
            tolerance,
            MinConjugateGradientTolerance,
            MaxConjugateGradientTolerance);
    }

    float World::GetConjugateGradientTolerance() const
    {
        return conjugateGradientSettings.tolerance;
    }

    void World::SetConjugateGradientMaxIterations(int maxIterations)
    {
        conjugateGradientSettings.maxIterations = std::clamp(
            maxIterations,
            MinConjugateGradientMaxIterations,
            MaxConjugateGradientMaxIterations);
    }

    int World::GetConjugateGradientMaxIterations() const
    {
        return conjugateGradientSettings.maxIterations;
    }

    int World::GetLastConjugateGradientIterations() const
    {
        return lastConjugateGradientResult.iterations;
    }

    float World::GetLastConjugateGradientResidualNorm() const
    {
        return lastConjugateGradientResult.residualNorm;
    }

    bool World::DidLastConjugateGradientSolveConverge() const
    {
        return lastConjugateGradientResult.converged;
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

    Vec3& World::GetNodePosition(int nodeIndex)
    {
        assert(nodeIndex >= 0 && nodeIndex < static_cast<int>(nodes.size()));
        const int stateIndex = nodes[static_cast<std::size_t>(nodeIndex)].stateIndex;
        assert(state.IsValidStateIndex(stateIndex));
        return state.positions[static_cast<std::size_t>(stateIndex)];
    }

    const Vec3& World::GetNodePosition(int nodeIndex) const
    {
        assert(nodeIndex >= 0 && nodeIndex < static_cast<int>(nodes.size()));
        const int stateIndex = nodes[static_cast<std::size_t>(nodeIndex)].stateIndex;
        assert(state.IsValidStateIndex(stateIndex));
        return state.positions[static_cast<std::size_t>(stateIndex)];
    }

    Vec3& World::GetNodeVelocity(int nodeIndex)
    {
        assert(nodeIndex >= 0 && nodeIndex < static_cast<int>(nodes.size()));
        const int stateIndex = nodes[static_cast<std::size_t>(nodeIndex)].stateIndex;
        assert(state.IsValidStateIndex(stateIndex));
        return state.velocities[static_cast<std::size_t>(stateIndex)];
    }

    const Vec3& World::GetNodeVelocity(int nodeIndex) const
    {
        assert(nodeIndex >= 0 && nodeIndex < static_cast<int>(nodes.size()));
        const int stateIndex = nodes[static_cast<std::size_t>(nodeIndex)].stateIndex;
        assert(state.IsValidStateIndex(stateIndex));
        return state.velocities[static_cast<std::size_t>(stateIndex)];
    }

    float& World::GetNodeMass(int nodeIndex)
    {
        assert(nodeIndex >= 0 && nodeIndex < static_cast<int>(nodes.size()));
        const int stateIndex = nodes[static_cast<std::size_t>(nodeIndex)].stateIndex;
        assert(state.IsValidStateIndex(stateIndex));
        return state.masses[static_cast<std::size_t>(stateIndex)];
    }

    float World::GetNodeMass(int nodeIndex) const
    {
        assert(nodeIndex >= 0 && nodeIndex < static_cast<int>(nodes.size()));
        const int stateIndex = nodes[static_cast<std::size_t>(nodeIndex)].stateIndex;
        assert(state.IsValidStateIndex(stateIndex));
        return state.masses[static_cast<std::size_t>(stateIndex)];
    }

    bool World::NodeHasRotation(int nodeIndex) const
    {
        assert(nodeIndex >= 0 && nodeIndex < static_cast<int>(nodes.size()));
        return nodes[static_cast<std::size_t>(nodeIndex)].hasRotation;
    }

    void World::SetNodeHasRotation(int nodeIndex, bool hasRotation)
    {
        assert(nodeIndex >= 0 && nodeIndex < static_cast<int>(nodes.size()));
        nodes[static_cast<std::size_t>(nodeIndex)].hasRotation = hasRotation;
    }

    Quaternion& World::GetNodeOrientation(int nodeIndex)
    {
        assert(nodeIndex >= 0 && nodeIndex < static_cast<int>(nodes.size()));
        const int stateIndex = nodes[static_cast<std::size_t>(nodeIndex)].stateIndex;
        assert(state.IsValidStateIndex(stateIndex));
        return state.orientations[static_cast<std::size_t>(stateIndex)];
    }

    const Quaternion& World::GetNodeOrientation(int nodeIndex) const
    {
        assert(nodeIndex >= 0 && nodeIndex < static_cast<int>(nodes.size()));
        const int stateIndex = nodes[static_cast<std::size_t>(nodeIndex)].stateIndex;
        assert(state.IsValidStateIndex(stateIndex));
        return state.orientations[static_cast<std::size_t>(stateIndex)];
    }

    void World::SetNodeOrientation(int nodeIndex, const Quaternion& orientation)
    {
        GetNodeOrientation(nodeIndex) = orientation;
    }

    Vec3& World::GetNodeAngularVelocity(int nodeIndex)
    {
        assert(nodeIndex >= 0 && nodeIndex < static_cast<int>(nodes.size()));
        const int stateIndex = nodes[static_cast<std::size_t>(nodeIndex)].stateIndex;
        assert(state.IsValidStateIndex(stateIndex));
        return state.angularVelocities[static_cast<std::size_t>(stateIndex)];
    }

    const Vec3& World::GetNodeAngularVelocity(int nodeIndex) const
    {
        assert(nodeIndex >= 0 && nodeIndex < static_cast<int>(nodes.size()));
        const int stateIndex = nodes[static_cast<std::size_t>(nodeIndex)].stateIndex;
        assert(state.IsValidStateIndex(stateIndex));
        return state.angularVelocities[static_cast<std::size_t>(stateIndex)];
    }

    void World::SetNodeAngularVelocity(int nodeIndex, const Vec3& angularVelocity)
    {
        GetNodeAngularVelocity(nodeIndex) = angularVelocity;
    }

    Vec3& World::GetNodeTorque(int nodeIndex)
    {
        assert(nodeIndex >= 0 && nodeIndex < static_cast<int>(nodes.size()));
        const int stateIndex = nodes[static_cast<std::size_t>(nodeIndex)].stateIndex;
        assert(state.IsValidStateIndex(stateIndex));
        return state.torques[static_cast<std::size_t>(stateIndex)];
    }

    const Vec3& World::GetNodeTorque(int nodeIndex) const
    {
        assert(nodeIndex >= 0 && nodeIndex < static_cast<int>(nodes.size()));
        const int stateIndex = nodes[static_cast<std::size_t>(nodeIndex)].stateIndex;
        assert(state.IsValidStateIndex(stateIndex));
        return state.torques[static_cast<std::size_t>(stateIndex)];
    }

    void World::SetNodeTorque(int nodeIndex, const Vec3& torque)
    {
        GetNodeTorque(nodeIndex) = torque;
    }

    Mat3& World::GetNodeInverseInertia(int nodeIndex)
    {
        assert(nodeIndex >= 0 && nodeIndex < static_cast<int>(nodes.size()));
        const int stateIndex = nodes[static_cast<std::size_t>(nodeIndex)].stateIndex;
        assert(state.IsValidStateIndex(stateIndex));
        return state.inverseInertias[static_cast<std::size_t>(stateIndex)];
    }

    const Mat3& World::GetNodeInverseInertia(int nodeIndex) const
    {
        assert(nodeIndex >= 0 && nodeIndex < static_cast<int>(nodes.size()));
        const int stateIndex = nodes[static_cast<std::size_t>(nodeIndex)].stateIndex;
        assert(state.IsValidStateIndex(stateIndex));
        return state.inverseInertias[static_cast<std::size_t>(stateIndex)];
    }

    void World::SetNodeInverseInertia(int nodeIndex, const Mat3& inverseInertia)
    {
        GetNodeInverseInertia(nodeIndex) = inverseInertia;
    }

    void World::SetNodePosition(int index, const Vec3& position)
    {
        assert(index >= 0 && index < static_cast<int>(nodes.size()));
        GetNodePosition(index) = position;
        GetNodeVelocity(index) = Vec3{};
        GetNodeAngularVelocity(index) = Vec3{};
        GetNodeTorque(index) = Vec3{};
    }

    void World::SetNodeFixed(int nodeIndex, bool fixed)
    {
        assert(nodeIndex >= 0 && nodeIndex < static_cast<int>(nodes.size()));
        Node& node = nodes[static_cast<std::size_t>(nodeIndex)];
        if (fixed)
        {
            node.fixed = true;
            GetNodeVelocity(nodeIndex) = Vec3{};
            GetNodeAngularVelocity(nodeIndex) = Vec3{};
            GetNodeTorque(nodeIndex) = Vec3{};
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

    ComponentHandle World::GetComponentHandleByIndex(int index) const
    {
        if (index < 0 || index >= static_cast<int>(components.size()) ||
            components[static_cast<std::size_t>(index)] == nullptr)
        {
            return ComponentHandle{};
        }

        return ComponentHandle{
            static_cast<std::uint32_t>(index),
            componentGenerations[static_cast<std::size_t>(index)]};
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
            Component* addedComponent = components[slotIndex].get();
            RegisterComponentForExecution(addedComponent);
            for (PhysicsEventType eventType : addedComponent->listenedEvents)
            {
                eventSystem.Subscribe(addedComponent, eventType);
            }
            return ComponentHandle{slotIndex, componentGenerations[slotIndex]};
        }

        components.push_back(std::move(component));
        componentGenerations.push_back(1u);
        Component* addedComponent = components.back().get();
        RegisterComponentForExecution(addedComponent);
        for (PhysicsEventType eventType : addedComponent->listenedEvents)
        {
            eventSystem.Subscribe(addedComponent, eventType);
        }

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

    void World::RegisterComponentForExecution(Component* component)
    {
        if (component == nullptr)
        {
            return;
        }

        orderedComponents.emplace(
            component->GetExecutionPriority(),
            component);
    }

    void World::UnregisterComponentFromExecution(Component* component)
    {
        if (component == nullptr)
        {
            return;
        }

        for (auto iterator = orderedComponents.begin();
             iterator != orderedComponents.end();)
        {
            if (iterator->second == component)
            {
                iterator = orderedComponents.erase(iterator);
            }
            else
            {
                ++iterator;
            }
        }
    }

    void World::PreUpdateComponents(float frameDt)
    {
        for (const auto& entry : orderedComponents)
        {
            Component* component = entry.second;
            if (component != nullptr && component->active)
            {
                component->PreUpdate(*this, frameDt);
            }
        }
    }

    void World::PostUpdateComponents(float frameDt)
    {
        for (const auto& entry : orderedComponents)
        {
            Component* component = entry.second;
            if (component != nullptr && component->active)
            {
                component->PostUpdate(*this, frameDt);
            }
        }
    }

    void World::BuildSolverData(SolverData& solverData, float dt)
    {
        AssembleComponentSystems(solverData, dt);
        solverData.AssembleMasses(static_cast<int>(nodes.size()));
        ValidateNodeMasses(solverData);
        AddGravityForces(solverData);
        AssembleConnectionSystems(solverData, dt);
    }

    void World::ValidateNodeMasses(const SolverData& solverData) const
    {
        for (int i = 0; i < static_cast<int>(nodes.size()); ++i)
        {
            const Node& node = nodes[static_cast<std::size_t>(i)];
            if (!node.active || node.fixed)
            {
                continue;
            }

            const float mass = solverData.GetAssembledMassForNode(i);
            if (std::isfinite(mass) && mass > 0.0f)
            {
                continue;
            }

            throw std::runtime_error(
                "Non-fixed node has no physical mass contribution. nodeIndex=" +
                std::to_string(i));
        }
    }

    void World::AddGravityForces(SolverData& solverData)
    {
        for (int i = 0; i < static_cast<int>(nodes.size()); ++i)
        {
            const Node& node = nodes[static_cast<std::size_t>(i)];
            if (!node.active || node.fixed)
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

    void World::AssembleComponentSystems(SolverData& solverData, float dt)
    {
        for (const auto& entry : orderedComponents)
        {
            Component* component = entry.second;
            if (component != nullptr && component->active)
            {
                component->UpdateSystem(*this, solverData, dt);
            }
        }
    }

    void World::PrecomputeSolve(SolverData& solverData, float dt)
    {
        solverData.ClearTransientState();
        BuildSolverData(solverData, dt);
        solverData.PrecomputeImplicitSolve(nodes, state.velocities, dt);
    }

    bool World::SolveImplicitLinearSystem(SolverData& solverData, float dt)
    {
        (void)dt;
        const bool solved =
            solverData.SolveImplicitLinearSystem(conjugateGradientSettings);
        const LinearSolveResult& result = solverData.GetLastLinearSolveResult();
        lastConjugateGradientResult.iterations = result.iterations;
        lastConjugateGradientResult.residualNorm = result.residualNorm;
        lastConjugateGradientResult.converged = result.converged;
        return solved;
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
            const Vec3& position = GetNodePosition(nodeIndex);
            const Vec3& velocity = GetNodeVelocity(nodeIndex);
            const std::vector<float>& deltaVelocity = solverData.GetDeltaVelocity();
            if (baseDof + 2 >= static_cast<int>(deltaVelocity.size()))
            {
                return false;
            }

            const Vec3 velocityChange{
                deltaVelocity[static_cast<std::size_t>(baseDof + 0)],
                deltaVelocity[static_cast<std::size_t>(baseDof + 1)],
                deltaVelocity[static_cast<std::size_t>(baseDof + 2)]};
            if (!IsFinite(position) || !IsFinite(velocity) || !IsFinite(velocityChange))
            {
                return false;
            }

            const Vec3 updatedVelocity = velocity + velocityChange;
            const Vec3 updatedPosition = position + updatedVelocity * dt;
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

            GetNodeVelocity(nodeIndex) = updatedVelocities[static_cast<std::size_t>(nodeIndex)];
            GetNodePosition(nodeIndex) = updatedPositions[static_cast<std::size_t>(nodeIndex)];
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
            if (!node.active ||
                node.fixed ||
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
            Vec3& velocity = GetNodeVelocity(i);
            Vec3& position = GetNodePosition(i);
            velocity += acceleration * dt;
            position += velocity * dt;
        }
    }

    void World::MarkSubstepConnectionBegin()
    {
        firstSubstepConnectionIndex = transientConnections.size();
    }

    void World::ClearSubstepConnections()
    {
        if (firstSubstepConnectionIndex > transientConnections.size())
        {
            firstSubstepConnectionIndex = transientConnections.size();
        }

        transientConnections.resize(firstSubstepConnectionIndex);
    }

    void World::ClearFrameConnections()
    {
        transientConnections.clear();
        firstSubstepConnectionIndex = 0u;
    }

    bool World::HasValidNodeIndices(const PointConnection& connection) const
    {
        const int nodeCount = static_cast<int>(nodes.size());
        return connection.node0 >= 0 && connection.node0 < nodeCount &&
            connection.node1 >= 0 && connection.node1 < nodeCount &&
            connection.node2 >= 0 && connection.node2 < nodeCount &&
            connection.node3 >= 0 && connection.node3 < nodeCount &&
            nodes[static_cast<std::size_t>(connection.node0)].active &&
            nodes[static_cast<std::size_t>(connection.node1)].active &&
            nodes[static_cast<std::size_t>(connection.node2)].active &&
            nodes[static_cast<std::size_t>(connection.node3)].active;
    }
}
