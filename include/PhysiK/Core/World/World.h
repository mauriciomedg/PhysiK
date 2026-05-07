#pragma once

#include <memory>
#include <vector>

#include "PhysiK/Components/Component.h"
#include "PhysiK/Components/TetMeshComponent.h"
#include "PhysiK/Math/Vec3.h"
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

        void AddPointConnection(const PointConnection& connection);

        void SetSubstepCount(int count);
        int GetSubstepCount() const;

        Node& GetNode(int index);
        const Node& GetNode(int index) const;

        const std::vector<PointConnection>& GetPointConnections() const;

    private:
        void ApplyPointConnectionForces();
        void Integrate(float dt);
        void ClearForces();
        bool HasValidNodeIndices(const PointConnection& connection) const;

        std::vector<Node> nodes;
        std::vector<Tet> tets;

        std::vector<std::unique_ptr<Component>> components;
        std::vector<TetMeshComponent*> tetMeshes;

        std::vector<PointConnection> pointConnections;

        int substepCount = 1;
    };
}
