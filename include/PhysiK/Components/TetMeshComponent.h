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
        std::vector<Vec3> currentNodePositions;
        std::vector<Tet> tets;

        bool topologyDirty = false;

        void SetGeometry(
            const Vec3* positions,
            int nodeCount,
            const int* tetLocalNodeIndices,
            int tetCount);

        int GetNodeCount() const;
        int GetTetCount() const;
        int GetTetNodeIndex(int tetIndex, int cornerIndex) const;
        const Vec3& GetLocalRestPosition(int localNodeIndex) const;
        const Vec3& GetLocalCurrentPosition(int localNodeIndex) const;
        virtual int GetWorldNodeIndex(int localNodeIndex) const;
        virtual void SetLocalCurrentPosition(int localNodeIndex, const Vec3& position);

        bool IsTetActive(int tetIndex) const;
        void SetTetActive(int tetIndex, bool active);
        void DeactivateTet(int tetIndex);
        int GetActiveTetCount() const;

        void PostUpdate(World& world, float dt) override;
    };
}
