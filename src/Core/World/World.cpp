#include "PhysiK/Core/World/World.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <memory>

#include "PhysiK/Core/Solvers/Linear/ConjugateGradientSolver.h"
#include "PhysiK/Math/SparseBlockMatrix.h"

namespace PhysiK
{
    struct World::ImplicitSolveData
    {
        std::vector<float> masses;
        std::vector<int> nodeToDynamicBlock;
        int dynamicBlockCount = 0;
        std::vector<float> rhs;
        SparseBlockMatrix matrix;
        std::vector<float> deltaVelocity;
    };

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

        bool IsFinite(const Mat3& matrix)
        {
            return IsFinite(matrix.columns[0]) &&
                IsFinite(matrix.columns[1]) &&
                IsFinite(matrix.columns[2]);
        }

        std::vector<float> BuildNodeMasses(const SolverData& solverData, int nodeCount)
        {
            std::vector<float> masses(static_cast<std::size_t>(nodeCount), 0.0f);
            for (const SolverData::NodeMass& nodeMass : solverData.GetNodeMasses())
            {
                if (nodeMass.node < 0 || nodeMass.node >= nodeCount ||
                    !std::isfinite(nodeMass.mass) || nodeMass.mass <= 0.0f)
                {
                    continue;
                }

                masses[static_cast<std::size_t>(nodeMass.node)] += nodeMass.mass;
            }

            return masses;
        }

        std::vector<bool> BuildNodeMassContributions(const SolverData& solverData, int nodeCount)
        {
            std::vector<bool> contributed(static_cast<std::size_t>(nodeCount), false);
            for (const SolverData::NodeMass& nodeMass : solverData.GetNodeMasses())
            {
                if (nodeMass.node < 0 || nodeMass.node >= nodeCount)
                {
                    continue;
                }

                contributed[static_cast<std::size_t>(nodeMass.node)] = true;
            }

            return contributed;
        }

        std::vector<Vec3> BuildNodeForces(const SolverData& solverData, int nodeCount)
        {
            std::vector<Vec3> forces(static_cast<std::size_t>(nodeCount));
            for (const SolverData::NodeForce& nodeForce : solverData.GetNodeForces())
            {
                if (nodeForce.node < 0 || nodeForce.node >= nodeCount ||
                    !IsFinite(nodeForce.force))
                {
                    continue;
                }

                forces[static_cast<std::size_t>(nodeForce.node)] += nodeForce.force;
            }

            return forces;
        }

        Mat3 ScaledIdentity(float value)
        {
            return Mat3::FromColumns(
                Vec3{value, 0.0f, 0.0f},
                Vec3{0.0f, value, 0.0f},
                Vec3{0.0f, 0.0f, value});
        }

        Mat3 Scale(const Mat3& matrix, float value)
        {
            return Mat3::FromColumns(
                matrix.columns[0] * value,
                matrix.columns[1] * value,
                matrix.columns[2] * value);
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
            BuildSolverData(solverData, substepDt);

            if (solverMode == SolverMode::ImplicitEuler)
            {
                SolveImplicitEuler(solverData, substepDt);
            }
            else
            {
                ApplyExplicitForces(solverData, substepDt);
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

        solverData.Clear();
        AssembleComponentSystems(solverData, dt);
        AddDefaultNodeMasses(solverData);
        AddGravityForces(solverData);
        AssembleConnectionSystems(solverData, dt);
    }

    void World::AddDefaultNodeMasses(SolverData& solverData)
    {
        const std::vector<float> masses =
            BuildNodeMasses(solverData, static_cast<int>(nodes.size()));
        const std::vector<bool> massContributed =
            BuildNodeMassContributions(solverData, static_cast<int>(nodes.size()));

        for (int i = 0; i < static_cast<int>(nodes.size()); ++i)
        {
            if (nodes[static_cast<std::size_t>(i)].fixed ||
                massContributed[static_cast<std::size_t>(i)] ||
                masses[static_cast<std::size_t>(i)] > 0.0f)
            {
                continue;
            }

            solverData.AddNodeMass(i, 1.0f);
        }
    }

    void World::AddGravityForces(SolverData& solverData)
    {
        const std::vector<float> masses =
            BuildNodeMasses(solverData, static_cast<int>(nodes.size()));

        for (int i = 0; i < static_cast<int>(nodes.size()); ++i)
        {
            if (nodes[static_cast<std::size_t>(i)].fixed)
            {
                continue;
            }

            const float mass = masses[static_cast<std::size_t>(i)];
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

    void World::ApplyExplicitForces(SolverData& solverData, float dt)
    {
        const std::vector<float> masses =
            BuildNodeMasses(solverData, static_cast<int>(nodes.size()));
        const std::vector<Vec3> forces =
            BuildNodeForces(solverData, static_cast<int>(nodes.size()));

        for (int i = 0; i < static_cast<int>(nodes.size()); ++i)
        {
            Node& node = nodes[static_cast<std::size_t>(i)];
            if (node.fixed)
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

    void World::SolveImplicitEuler(SolverData& solverData, float dt)
    {
        ImplicitSolveData solveData;
        if (!PrecomputeSolve(solverData, dt, solveData))
        {
            return;
        }

        if (!SolveImplicitLinearSystem(solveData))
        {
            return;
        }

        if (!IntegrateImplicitEuler(solveData, dt))
        {
            return;
        }
    }

    bool World::PrecomputeSolve(
        const SolverData& solverData,
        float dt,
        ImplicitSolveData& solveData) const
    {
        if (dt <= 0.0f)
        {
            return false;
        }

        solveData = ImplicitSolveData{};
        solveData.masses = BuildNodeMasses(solverData, static_cast<int>(nodes.size()));
        solveData.nodeToDynamicBlock.assign(nodes.size(), -1);

        for (int nodeIndex = 0; nodeIndex < static_cast<int>(nodes.size()); ++nodeIndex)
        {
            const float mass = solveData.masses[static_cast<std::size_t>(nodeIndex)];
            if (!nodes[static_cast<std::size_t>(nodeIndex)].fixed &&
                std::isfinite(mass) && mass > 0.0f)
            {
                solveData.nodeToDynamicBlock[static_cast<std::size_t>(nodeIndex)] =
                    solveData.dynamicBlockCount;
                ++solveData.dynamicBlockCount;
            }
        }

        if (solveData.dynamicBlockCount == 0)
        {
            return false;
        }

        const std::size_t dimension =
            static_cast<std::size_t>(solveData.dynamicBlockCount * 3);
        solveData.rhs.assign(dimension, 0.0f);

        std::vector<std::pair<int, int>> blockCoordinates;
        blockCoordinates.reserve(
            static_cast<std::size_t>(solveData.dynamicBlockCount) +
            solverData.GetStiffnessBlocks().size());

        for (int nodeIndex = 0; nodeIndex < static_cast<int>(nodes.size()); ++nodeIndex)
        {
            const int dynamicBlock =
                solveData.nodeToDynamicBlock[static_cast<std::size_t>(nodeIndex)];
            if (dynamicBlock < 0)
            {
                continue;
            }

            const float mass = solveData.masses[static_cast<std::size_t>(nodeIndex)];
            if (!IsFinite(mass))
            {
                return false;
            }

            blockCoordinates.push_back({dynamicBlock, dynamicBlock});
        }

        const float stiffnessScale = dt * dt;
        for (const SolverData::StiffnessBlock& block : solverData.GetStiffnessBlocks())
        {
            if (block.nodeA < 0 || block.nodeA >= static_cast<int>(nodes.size()) ||
                block.nodeB < 0 || block.nodeB >= static_cast<int>(nodes.size()))
            {
                continue;
            }

            const int rowBlock =
                solveData.nodeToDynamicBlock[static_cast<std::size_t>(block.nodeA)];
            if (rowBlock < 0)
            {
                continue;
            }

            if (!IsFinite(block.block))
            {
                return false;
            }

            const Vec3& columnVelocity = nodes[static_cast<std::size_t>(block.nodeB)].velocity;
            if (!IsFinite(columnVelocity))
            {
                return false;
            }

            const Vec3 stiffnessVelocity = block.block * columnVelocity;
            if (!IsFinite(stiffnessVelocity))
            {
                return false;
            }

            const int rowBase = rowBlock * 3;
            solveData.rhs[static_cast<std::size_t>(rowBase + 0)] -=
                stiffnessScale * stiffnessVelocity.x;
            solveData.rhs[static_cast<std::size_t>(rowBase + 1)] -=
                stiffnessScale * stiffnessVelocity.y;
            solveData.rhs[static_cast<std::size_t>(rowBase + 2)] -=
                stiffnessScale * stiffnessVelocity.z;

            const int columnBlock =
                solveData.nodeToDynamicBlock[static_cast<std::size_t>(block.nodeB)];
            if (columnBlock < 0)
            {
                continue;
            }

            blockCoordinates.push_back({rowBlock, columnBlock});
        }

        solveData.matrix.BuildPattern(solveData.dynamicBlockCount, blockCoordinates);

        for (int nodeIndex = 0; nodeIndex < static_cast<int>(nodes.size()); ++nodeIndex)
        {
            const int dynamicBlock =
                solveData.nodeToDynamicBlock[static_cast<std::size_t>(nodeIndex)];
            if (dynamicBlock < 0)
            {
                continue;
            }

            const float mass = solveData.masses[static_cast<std::size_t>(nodeIndex)];
            if (!solveData.matrix.AddBlock(dynamicBlock, dynamicBlock, ScaledIdentity(mass)))
            {
                return false;
            }
        }

        for (const SolverData::StiffnessBlock& block : solverData.GetStiffnessBlocks())
        {
            if (block.nodeA < 0 || block.nodeA >= static_cast<int>(nodes.size()) ||
                block.nodeB < 0 || block.nodeB >= static_cast<int>(nodes.size()) ||
                !IsFinite(block.block))
            {
                continue;
            }

            const int rowBlock =
                solveData.nodeToDynamicBlock[static_cast<std::size_t>(block.nodeA)];
            const int columnBlock =
                solveData.nodeToDynamicBlock[static_cast<std::size_t>(block.nodeB)];
            if (rowBlock < 0 || columnBlock < 0)
            {
                continue;
            }

            if (!solveData.matrix.AddBlock(rowBlock, columnBlock, Scale(block.block, stiffnessScale)))
            {
                return false;
            }
        }

        for (const SolverData::NodeForce& nodeForce : solverData.GetNodeForces())
        {
            if (nodeForce.node < 0 || nodeForce.node >= static_cast<int>(nodes.size()))
            {
                continue;
            }

            if (!IsFinite(nodeForce.force))
            {
                return false;
            }

            const int dynamicBlock =
                solveData.nodeToDynamicBlock[static_cast<std::size_t>(nodeForce.node)];
            if (dynamicBlock < 0)
            {
                continue;
            }

            const int baseDof = dynamicBlock * 3;
            solveData.rhs[static_cast<std::size_t>(baseDof + 0)] += dt * nodeForce.force.x;
            solveData.rhs[static_cast<std::size_t>(baseDof + 1)] += dt * nodeForce.force.y;
            solveData.rhs[static_cast<std::size_t>(baseDof + 2)] += dt * nodeForce.force.z;
        }

        return true;
    }

    bool World::SolveImplicitLinearSystem(ImplicitSolveData& solveData) const
    {
        const std::size_t dimension =
            static_cast<std::size_t>(solveData.dynamicBlockCount * 3);
        if (dimension == 0)
        {
            return false;
        }

        ConjugateGradientSettings cgSettings;
        cgSettings.maxIterations = std::max(128, solveData.dynamicBlockCount * 12);
        cgSettings.tolerance = 1.0e-5f;
        cgSettings.useJacobiPreconditioner = true;

        const ConjugateGradientResult cgResult =
            SolveConjugateGradient(
                solveData.matrix,
                solveData.rhs,
                solveData.deltaVelocity,
                cgSettings);
        return cgResult.converged && solveData.deltaVelocity.size() == dimension;
    }

    bool World::IntegrateImplicitEuler(const ImplicitSolveData& solveData, float dt)
    {
        std::vector<Vec3> updatedVelocities(nodes.size());
        std::vector<Vec3> updatedPositions(nodes.size());

        for (int nodeIndex = 0; nodeIndex < static_cast<int>(nodes.size()); ++nodeIndex)
        {
            const int dynamicBlock =
                solveData.nodeToDynamicBlock[static_cast<std::size_t>(nodeIndex)];
            if (dynamicBlock < 0)
            {
                continue;
            }

            const int baseDof = dynamicBlock * 3;
            const Node& node = nodes[static_cast<std::size_t>(nodeIndex)];
            const Vec3 velocityChange{
                solveData.deltaVelocity[static_cast<std::size_t>(baseDof + 0)],
                solveData.deltaVelocity[static_cast<std::size_t>(baseDof + 1)],
                solveData.deltaVelocity[static_cast<std::size_t>(baseDof + 2)]};
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
            if (solveData.nodeToDynamicBlock[static_cast<std::size_t>(nodeIndex)] < 0)
            {
                continue;
            }

            Node& node = nodes[static_cast<std::size_t>(nodeIndex)];
            node.velocity = updatedVelocities[static_cast<std::size_t>(nodeIndex)];
            node.position = updatedPositions[static_cast<std::size_t>(nodeIndex)];
        }

        return true;
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
