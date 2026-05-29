#pragma once

#include <memory>
#include <vector>

#include "PhysiK/API/PhysiKAPI.h"
#include "PhysiK/Components/Component.h"
#include "PhysiK/Math/Vec3.h"
#include "PhysiK/PhysicsData/Tet.h"

namespace PhysiK
{
    class PHYSIK_API TetMeshComponent : public Component
    {
    public:
        std::vector<Vec3> restNodePositions;
        std::vector<Vec3> nodePositions;
        std::vector<Tet> tets;

        bool topologyDirty = false;

        void SetGeometry(
            const Vec3* positions,
            int nodeCount,
            const int* tetLocalNodeIndices,
            int tetCount);
        void SetGeometry(
            const Vec3* positions,
            int nodeCount,
            const int* tetLocalNodeIndices,
            int tetCount,
            std::vector<int>* outOldNodeToNewNode,
            std::vector<int>* outNewNodeToFirstOldNode);

        int GetNodeCount() const;
        int GetTetCount() const;
        int GetTetNodeIndex(int tetIndex, int cornerIndex) const;
        const Vec3& GetLocalRestPosition(int localNodeIndex) const;
        const Vec3& GetLocalCurrentPosition(int localNodeIndex) const;
        virtual int GetGlobalNodeIndex(int localNodeIndex) const;
        virtual void SetLocalCurrentPosition(int localNodeIndex, const Vec3& position);

        bool IsTetActive(int tetIndex) const;
        virtual bool SetTetActive(int tetIndex, bool active);
        virtual bool DeactivateTet(int tetIndex);
        int GetActiveTetCount() const;

        void PostUpdate(World& world, float dt) override;
    };
}
