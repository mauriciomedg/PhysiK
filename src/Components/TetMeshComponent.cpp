#include "PhysiK/Components/TetMeshComponent.h"

#include <cstddef>

#include "PhysiK/Core/World/World.h"

namespace PhysiK
{
    namespace
    {
        const Vec3 ZeroVec3;

        Tet MakeTet(int node0, int node1, int node2, int node3)
        {
            Tet tet;
            tet.node0 = node0;
            tet.node1 = node1;
            tet.node2 = node2;
            tet.node3 = node3;
            return tet;
        }
    }

    void TetMeshComponent::SetGeometry(
        const Vec3* positions,
        int nodeCount,
        const int* tetLocalNodeIndices,
        int tetCount)
    {
        restNodePositions.clear();
        nodePositions.clear();
        tets.clear();

        if (positions != nullptr && nodeCount > 0)
        {
            restNodePositions.assign(positions, positions + nodeCount);
            nodePositions = restNodePositions;
        }

        if (tetLocalNodeIndices != nullptr && tetCount > 0)
        {
            tets.reserve(static_cast<std::size_t>(tetCount));
            for (int tetIndex = 0; tetIndex < tetCount; ++tetIndex)
            {
                const int local0 = tetLocalNodeIndices[tetIndex * 4 + 0];
                const int local1 = tetLocalNodeIndices[tetIndex * 4 + 1];
                const int local2 = tetLocalNodeIndices[tetIndex * 4 + 2];
                const int local3 = tetLocalNodeIndices[tetIndex * 4 + 3];
                const int localNodeCount = static_cast<int>(restNodePositions.size());

                if (local0 < 0 || local0 >= localNodeCount ||
                    local1 < 0 || local1 >= localNodeCount ||
                    local2 < 0 || local2 >= localNodeCount ||
                    local3 < 0 || local3 >= localNodeCount)
                {
                    continue;
                }

                tets.push_back(MakeTet(local0, local1, local2, local3));
            }
        }
    }

    int TetMeshComponent::GetNodeCount() const
    {
        return static_cast<int>(restNodePositions.size());
    }

    int TetMeshComponent::GetTetCount() const
    {
        return static_cast<int>(tets.size());
    }

    int TetMeshComponent::GetTetNodeIndex(int tetIndex, int cornerIndex) const
    {
        if (tetIndex < 0 || tetIndex >= static_cast<int>(tets.size()))
        {
            return -1;
        }

        const Tet& tet = tets[static_cast<std::size_t>(tetIndex)];
        switch (cornerIndex)
        {
        case 0:
            return tet.node0;
        case 1:
            return tet.node1;
        case 2:
            return tet.node2;
        case 3:
            return tet.node3;
        default:
            return -1;
        }
    }

    const Vec3& TetMeshComponent::GetLocalRestPosition(int localNodeIndex) const
    {
        if (localNodeIndex < 0 ||
            localNodeIndex >= static_cast<int>(restNodePositions.size()))
        {
            return ZeroVec3;
        }

        return restNodePositions[static_cast<std::size_t>(localNodeIndex)];
    }

    const Vec3& TetMeshComponent::GetLocalCurrentPosition(int localNodeIndex) const
    {
        if (localNodeIndex < 0 ||
            localNodeIndex >= static_cast<int>(nodePositions.size()))
        {
            return ZeroVec3;
        }

        return nodePositions[static_cast<std::size_t>(localNodeIndex)];
    }

    int TetMeshComponent::GetGlobalNodeIndex(int localNodeIndex) const
    {
        (void)localNodeIndex;
        return -1;
    }

    void TetMeshComponent::SetLocalCurrentPosition(
        int localNodeIndex,
        const Vec3& position)
    {
        if (localNodeIndex < 0 ||
            localNodeIndex >= static_cast<int>(nodePositions.size()))
        {
            return;
        }

        nodePositions[static_cast<std::size_t>(localNodeIndex)] = position;
    }

    bool TetMeshComponent::IsTetActive(int tetIndex) const
    {
        if (tetIndex < 0 || tetIndex >= static_cast<int>(tets.size()))
        {
            return false;
        }

        return tets[static_cast<std::size_t>(tetIndex)].active;
    }

    void TetMeshComponent::SetTetActive(int tetIndex, bool active)
    {
        if (tetIndex < 0 || tetIndex >= static_cast<int>(tets.size()))
        {
            return;
        }

        tets[static_cast<std::size_t>(tetIndex)].active = active;
        topologyDirty = true;
    }

    void TetMeshComponent::DeactivateTet(int tetIndex)
    {
        SetTetActive(tetIndex, false);
    }

    int TetMeshComponent::GetActiveTetCount() const
    {
        int activeCount = 0;
        for (const Tet& tet : tets)
        {
            if (tet.active)
            {
                ++activeCount;
            }
        }

        return activeCount;
    }

    void TetMeshComponent::PostUpdate(World& world, float dt)
    {
        (void)dt;

        if (!topologyDirty)
        {
            return;
        }

        PhysicsEvent event;
        event.type = PhysicsEventType::TetMeshTopologyChanged;
        event.world = &world;
        event.sender = this;
        world.EmitEvent(event);

        topologyDirty = false;
    }
}
