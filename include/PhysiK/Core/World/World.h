#pragma once

#include <memory>
#include <vector>

#include "PhysiK/API/Handles.h"
#include "PhysiK/Components/CollisionComponent.h"
#include "PhysiK/Components/CollisionSphereComponent.h"
#include "PhysiK/Components/Component.h"
#include "PhysiK/Components/TetMeshComponent.h"
#include "PhysiK/Core/Collision/CollisionDetectionEngine.h"
#include "PhysiK/Core/PhysicsConnections/LineConnection.h"
#include "PhysiK/Core/PhysicsConnections/PointConnection.h"
#include "PhysiK/Core/PhysicsConnections/RigidBodyConnection.h"
#include "PhysiK/Core/PhysicsConnections/RigidBodyOrientationConnection.h"
#include "PhysiK/Core/PhysicsConnections/SurfaceConnection.h"
#include "PhysiK/Core/Solvers/SolverData.h"
#include "PhysiK/Math/Vec3.h"
#include "PhysiK/PhysicsData/Contact.h"
#include "PhysiK/PhysicsData/Node.h"

namespace PhysiK
{
    class World
    {
    public:
        World();

        void Step(float frameDt);

        int AddNode(const Vec3& position, float inverseMass = 1.0f);
        int AddTet(int node0, int node1, int node2, int node3);

        ComponentHandle CreateTetMeshComponent(
            const int* nodeIndices,
            int nodeCount,
            const int* tetNodeIndices,
            int tetCount);
        ComponentHandle CreateCollisionSphereComponent(
            const Vec3& position,
            float radius);
        Component* GetComponent(ComponentHandle handle);
        const Component* GetComponent(ComponentHandle handle) const;
        void DestroyComponent(ComponentHandle handle);

        void AddPointConnection(const PointConnection& connection);
        void SetExternalLogicCallback(ExternalLogicCallback callback, void* userData);
        void ClearExternalLogicCallback();

        void SetSubstepCount(int count);
        int GetSubstepCount() const;
        void SetGravity(const Vec3& value);
        const Vec3& GetGravity() const;

        Node& GetNode(int index);
        const Node& GetNode(int index) const;
        const std::vector<Node>& GetNodes() const;
        void SetNodePosition(int index, const Vec3& position);
        const std::vector<TetMeshComponent*>& GetTetMeshes() const;

        const std::vector<PointConnection>& GetPointConnections() const;
        bool HasValidNodeIndices(const PointConnection& connection) const;

    private:
        struct ComponentSlot
        {
            std::unique_ptr<Component> component;
            std::uint32_t generation = 1u;
        };

        ComponentHandle StoreComponent(std::unique_ptr<Component> component);
        bool IsComponentHandleValid(ComponentHandle handle) const;
        void RemoveTypedComponentReferences(Component* component);
        void RunExternalLogic();
        void UpdateKinematicTargets();
        void AccumulateForces(float dt);
        void AddGravityForces(SolverData& solverData);
        void AddConnectionForces(SolverData& solverData, float dt);
        void AddCollisionForces(SolverData& solverData, float dt);
        void AddPhysicsModelForces(SolverData& solverData, float dt);
        void AddPointConnectionFromContact(const Contact& contact, SolverData& solverData, float dt);
        void Solve(SolverData& solverData, float dt);
        void Integrate(float dt);
        void ClearForces();
        void ClearTransientConnections();

        std::vector<Node> nodes;

        std::vector<ComponentSlot> componentSlots;
        std::vector<std::uint32_t> freeComponentSlots;
        std::vector<Component*> components;
        std::vector<TetMeshComponent*> tetMeshes;
        std::vector<CollisionComponent*> collisionComponents;

        std::vector<PointConnection> pointConnections;
        std::vector<SurfaceConnection> surfaceConnections;
        std::vector<LineConnection> lineConnections;
        std::vector<RigidBodyConnection> rigidBodyConnections;
        std::vector<RigidBodyOrientationConnection> rigidBodyOrientationConnections;

        CollisionDetectionEngine collisionDetectionEngine;

        ExternalLogicCallback externalLogicCallback = nullptr;
        void* externalLogicUserData = nullptr;

        Vec3 gravity;
        int substepCount = 1;
    };
}
