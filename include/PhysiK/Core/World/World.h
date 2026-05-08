#pragma once

#include <memory>
#include <vector>

#include "PhysiK/API/Handles.h"
#include "PhysiK/Components/CollisionComponent.h"
#include "PhysiK/Components/CollisionSphereComponent.h"
#include "PhysiK/Components/Component.h"
#include "PhysiK/Components/TetMeshComponent.h"
#include "PhysiK/Core/Collision/CollisionDetectionEngine.h"
#include "PhysiK/Core/Physics/FEM/FEMModel.h"
#include "PhysiK/Math/Vec3.h"
#include "PhysiK/PhysicsData/Contact.h"
#include "PhysiK/PhysicsData/Node.h"
#include "PhysiK/PhysicsData/PointConnection.h"
#include "PhysiK/PhysicsData/Tet.h"

namespace PhysiK
{
    class World
    {
    public:
        void Step(float frameDt);

        int AddNode(const Vec3& position, float inverseMass = 1.0f);
        int AddTet(int node0, int node1, int node2, int node3);

        ComponentHandle CreateTetMeshComponent(
            const int* nodeIndices,
            int nodeCount,
            const int* tetIndices,
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
        void SetNodePosition(int index, const Vec3& position);
        const std::vector<Tet>& GetTets() const;

        const std::vector<PointConnection>& GetPointConnections() const;

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
        void AccumulateForces();
        void AddGravityForces();
        void AddConnectionForces();
        void AddCollisionForces();
        void ApplyFEMForces();
        void AddPointConnectionFromContact(const Contact& contact);
        void AddPointConnectionForce(const PointConnection& connection);
        void Integrate(float dt);
        void ClearForces();
        bool HasValidNodeIndices(const PointConnection& connection) const;

        std::vector<Node> nodes;
        std::vector<Tet> tets;

        std::vector<ComponentSlot> componentSlots;
        std::vector<std::uint32_t> freeComponentSlots;
        std::vector<TetMeshComponent*> tetMeshes;
        std::vector<CollisionComponent*> collisionComponents;

        std::vector<PointConnection> pointConnections;

        CollisionDetectionEngine collisionDetectionEngine;

        ExternalLogicCallback externalLogicCallback = nullptr;
        void* externalLogicUserData = nullptr;

        Vec3 gravity;
        int substepCount = 1;
    };
}
