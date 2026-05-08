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

    World::World()
    {
        physicsModels.push_back(&femModel);
    }

    void World::Step(float frameDt)
    {
        if (frameDt <= 0.0f)
        {
            return;
        }

        RunExternalLogic();
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

    int World::AddTet(int node0, int node1, int node2, int node3)
    {
        tets.push_back(Tet{node0, node1, node2, node3});
        FEMModel::InitializeTetRestData(tets.back(), nodes);
        return static_cast<int>(tets.size()) - 1;
    }

    ComponentHandle World::CreateTetMeshComponent(
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

            for (int tetIndex : component->tetIndices)
            {
                if (tetIndex >= 0 && tetIndex < static_cast<int>(tets.size()))
                {
                    Tet& tet = tets[static_cast<std::size_t>(tetIndex)];
                    tet.stiffness = component->material.stiffness;
                    tet.damping = component->material.damping;
                    FEMModel::InitializeTetRestData(tet, nodes);
                }
            }
        }

        TetMeshComponent& componentRef = *component;
        tetMeshes.push_back(&componentRef);
        return StoreComponent(std::move(component));
    }

    ComponentHandle World::CreateCollisionSphereComponent(
        const Vec3& position,
        float radius)
    {
        auto component = std::make_unique<CollisionSphereComponent>();
        component->transform.position = position;
        component->radius = std::max(0.0f, radius);

        CollisionSphereComponent& componentRef = *component;
        collisionComponents.push_back(&componentRef);
        return StoreComponent(std::move(component));
    }

    Component* World::GetComponent(ComponentHandle handle)
    {
        if (!IsComponentHandleValid(handle))
        {
            return nullptr;
        }

        return componentSlots[handle.index].component.get();
    }

    const Component* World::GetComponent(ComponentHandle handle) const
    {
        if (!IsComponentHandleValid(handle))
        {
            return nullptr;
        }

        return componentSlots[handle.index].component.get();
    }

    void World::DestroyComponent(ComponentHandle handle)
    {
        if (!IsComponentHandleValid(handle))
        {
            return;
        }

        ComponentSlot& slot = componentSlots[handle.index];
        RemoveTypedComponentReferences(slot.component.get());
        slot.component.reset();
        ++slot.generation;

        if (slot.generation == 0u)
        {
            slot.generation = 1u;
        }

        freeComponentSlots.push_back(handle.index);
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

    const std::vector<Tet>& World::GetTets() const
    {
        return tets;
    }

    const std::vector<PointConnection>& World::GetPointConnections() const
    {
        return pointConnections;
    }

    ComponentHandle World::StoreComponent(std::unique_ptr<Component> component)
    {
        if (!freeComponentSlots.empty())
        {
            const std::uint32_t slotIndex = freeComponentSlots.back();
            freeComponentSlots.pop_back();
            ComponentSlot& slot = componentSlots[slotIndex];
            slot.component = std::move(component);
            return ComponentHandle{slotIndex, slot.generation};
        }

        ComponentSlot slot;
        slot.component = std::move(component);
        componentSlots.push_back(std::move(slot));

        return ComponentHandle{
            static_cast<std::uint32_t>(componentSlots.size() - 1u),
            componentSlots.back().generation};
    }

    bool World::IsComponentHandleValid(ComponentHandle handle) const
    {
        if (!handle.IsValid() || handle.index >= componentSlots.size())
        {
            return false;
        }

        const ComponentSlot& slot = componentSlots[handle.index];
        return slot.generation == handle.generation && slot.component != nullptr;
    }

    void World::RemoveTypedComponentReferences(Component* component)
    {
        if (component == nullptr)
        {
            return;
        }

        tetMeshes.erase(
            std::remove_if(
                tetMeshes.begin(),
                tetMeshes.end(),
                [component](const TetMeshComponent* tetMesh)
                {
                    return tetMesh == component;
                }),
            tetMeshes.end());

        if (auto* collision = dynamic_cast<CollisionComponent*>(component))
        {
            collisionComponents.erase(
                std::remove(collisionComponents.begin(), collisionComponents.end(), collision),
                collisionComponents.end());
        }
    }

    void World::RunExternalLogic()
    {
        if (externalLogicCallback != nullptr)
        {
            externalLogicCallback(static_cast<WorldHandle>(this), externalLogicUserData);
        }
    }

    void World::UpdateKinematicTargets()
    {
        for (CollisionComponent* component : collisionComponents)
        {
            if (component == nullptr)
            {
                continue;
            }

            Transform target;
            if (component->ConsumeKinematicTarget(target))
            {
                component->transform = target;
            }
        }
    }

    void World::AccumulateForces(float dt)
    {
        ClearForces();
        SolverData solverData;
        solverData.Clear();
        AddGravityForces(solverData);
        AddConnectionForces(solverData, dt);
        AddCollisionForces(solverData, dt);
        AddPhysicsModelForces(solverData, dt);
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
        for (PointConnection& connection : pointConnections)
        {
            connection.UpdateSystem(*this, solverData, dt);
        }

        for (SurfaceConnection& connection : surfaceConnections)
        {
            connection.UpdateSystem(*this, solverData, dt);
        }

        for (LineConnection& connection : lineConnections)
        {
            connection.UpdateSystem(*this, solverData, dt);
        }

        for (RigidBodyConnection& connection : rigidBodyConnections)
        {
            connection.UpdateSystem(*this, solverData, dt);
        }

        for (RigidBodyOrientationConnection& connection : rigidBodyOrientationConnections)
        {
            connection.UpdateSystem(*this, solverData, dt);
        }
    }

    void World::AddCollisionForces(SolverData& solverData, float dt)
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
                AddPointConnectionFromContact(contact, solverData, dt);
            }
        }
    }

    void World::AddPhysicsModelForces(SolverData& solverData, float dt)
    {
        for (PhysicsModel* model : physicsModels)
        {
            if (model != nullptr)
            {
                model->UpdateSystem(*this, solverData, dt);
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
        connection.node0 = contact.tetNode0;
        connection.node1 = contact.tetNode1;
        connection.node2 = contact.tetNode2;
        connection.node3 = contact.tetNode3;
        connection.barycentric = contact.barycentric;
        connection.targetPosition = contact.worldPoint + contact.normal * contact.penetrationDepth;
        connection.stiffness = contact.stiffness;
        connection.damping = contact.damping;
        AddPointConnection(connection);
        connection.UpdateSystem(*this, solverData, dt);
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
        pointConnections.clear();
        surfaceConnections.clear();
        lineConnections.clear();
        rigidBodyConnections.clear();
        rigidBodyOrientationConnections.clear();
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
