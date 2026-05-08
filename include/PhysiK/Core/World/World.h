#pragma once

#include <memory>
#include <vector>

#include "PhysiK/API/Handles.h"
#include "PhysiK/Components/CollisionComponent.h"
#include "PhysiK/Components/CollisionSphereComponent.h"
#include "PhysiK/Components/Component.h"
#include "PhysiK/Components/TetMeshComponent.h"
#include "PhysiK/Core/Collision/CollisionDetectionEngine.h"
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

        TetMeshComponent& CreateTetMeshComponent(
            const int* nodeIndices,
            int nodeCount,
            const int* tetIndices,
            int tetCount);
        CollisionSphereComponent& CreateCollisionSphereComponent(
            const Vec3& position,
            float radius);

        void AddPointConnection(const PointConnection& connection);
        void SetExternalLogicCallback(ExternalLogicCallback callback, void* userData);
        void ClearExternalLogicCallback();

        void SetSubstepCount(int count);
        int GetSubstepCount() const;

        Node& GetNode(int index);
        const Node& GetNode(int index) const;
        const std::vector<Tet>& GetTets() const;

        const std::vector<PointConnection>& GetPointConnections() const;
        bool IsCollisionComponent(const CollisionComponent* component) const;

    private:
        void RunExternalLogic();
        void UpdateKinematicTargets();
        void GenerateCollisionConnections();
        void AddPointConnectionFromContact(const Contact& contact);
        void ApplyPointConnectionForces();
        void Integrate(float dt);
        void ClearForces();
        bool HasValidNodeIndices(const PointConnection& connection) const;

        std::vector<Node> nodes;
        std::vector<Tet> tets;

        std::vector<std::unique_ptr<Component>> components;
        std::vector<TetMeshComponent*> tetMeshes;
        std::vector<CollisionComponent*> collisionComponents;

        std::vector<PointConnection> pointConnections;

        CollisionDetectionEngine collisionDetectionEngine;

        ExternalLogicCallback externalLogicCallback = nullptr;
        void* externalLogicUserData = nullptr;

        int substepCount = 1;
    };
}
