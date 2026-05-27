#include "PhysiK/Components/TopologyMeshComponent.h"

#include <cstddef>
#include <queue>

#include "PhysiK/Components/TetMeshComponent.h"
#include "PhysiK/Core/World/World.h"

namespace PhysiK
{
    namespace
    {
        int CountSharedNodes(
            const TetMeshComponent& tetMesh,
            int tetA,
            int tetB)
        {
            const int aNodes[4] = {
                tetMesh.GetTetNodeIndex(tetA, 0),
                tetMesh.GetTetNodeIndex(tetA, 1),
                tetMesh.GetTetNodeIndex(tetA, 2),
                tetMesh.GetTetNodeIndex(tetA, 3)};
            const int bNodes[4] = {
                tetMesh.GetTetNodeIndex(tetB, 0),
                tetMesh.GetTetNodeIndex(tetB, 1),
                tetMesh.GetTetNodeIndex(tetB, 2),
                tetMesh.GetTetNodeIndex(tetB, 3)};

            int sharedCount = 0;
            for (int aNode : aNodes)
            {
                for (int bNode : bNodes)
                {
                    if (aNode == bNode)
                    {
                        ++sharedCount;
                        break;
                    }
                }
            }

            return sharedCount;
        }

        bool AreTetNeighbors(
            const TetMeshComponent& tetMesh,
            int tetA,
            int tetB)
        {
            return CountSharedNodes(tetMesh, tetA, tetB) >= 3;
        }
    }

    TopologyMeshComponent::TopologyMeshComponent()
    {
        listenedEvents.push_back(PhysicsEventType::TetMeshTopologyChanged);
        emittedEvents.push_back(PhysicsEventType::TopologyMeshUpdated);
    }

    TopologyMeshComponent::TopologyMeshComponent(ComponentHandle hostTetMeshHandle)
        : TopologyMeshComponent()
    {
        this->hostTetMeshHandle = hostTetMeshHandle;
    }

    int TopologyMeshComponent::GetTetIslandId(int tetIndex) const
    {
        if (tetIndex < 0 || tetIndex >= static_cast<int>(tetIslandIds.size()))
        {
            return -1;
        }

        return tetIslandIds[static_cast<std::size_t>(tetIndex)];
    }

    int TopologyMeshComponent::GetIslandCount() const
    {
        return islandCount;
    }

    void TopologyMeshComponent::OnPhysicsEvent(const PhysicsEvent& event)
    {
        if (event.type != PhysicsEventType::TetMeshTopologyChanged ||
            event.world == nullptr)
        {
            return;
        }

        if (event.world->GetComponent(hostTetMeshHandle) != event.sender)
        {
            return;
        }

        topologyDirty = true;
    }

    void TopologyMeshComponent::PostUpdate(World& world, float dt)
    {
        (void)dt;

        if (!topologyDirty)
        {
            return;
        }

        RebuildTopology(world);
        topologyDirty = false;

        PhysicsEvent event;
        event.type = PhysicsEventType::TopologyMeshUpdated;
        event.world = &world;
        event.sender = this;
        world.EmitEvent(event);
    }

    void TopologyMeshComponent::RebuildTopology(World& world)
    {
        tetIslandIds.clear();
        islandCount = 0;

        const Component* hostComponent = world.GetComponent(hostTetMeshHandle);
        const TetMeshComponent* hostTetMesh =
            dynamic_cast<const TetMeshComponent*>(hostComponent);
        if (hostTetMesh == nullptr)
        {
            return;
        }

        const int tetCount = hostTetMesh->GetTetCount();
        tetIslandIds.assign(static_cast<std::size_t>(tetCount), -1);

        std::queue<int> pendingTets;
        for (int startTet = 0; startTet < tetCount; ++startTet)
        {
            if (!hostTetMesh->IsTetActive(startTet) ||
                tetIslandIds[static_cast<std::size_t>(startTet)] != -1)
            {
                continue;
            }

            const int islandId = islandCount;
            ++islandCount;
            tetIslandIds[static_cast<std::size_t>(startTet)] = islandId;
            pendingTets.push(startTet);

            while (!pendingTets.empty())
            {
                const int currentTet = pendingTets.front();
                pendingTets.pop();

                for (int candidateTet = 0;
                     candidateTet < tetCount;
                     ++candidateTet)
                {
                    if (!hostTetMesh->IsTetActive(candidateTet) ||
                        tetIslandIds[static_cast<std::size_t>(candidateTet)] != -1)
                    {
                        continue;
                    }

                    if (!AreTetNeighbors(
                            *hostTetMesh,
                            currentTet,
                            candidateTet))
                    {
                        continue;
                    }

                    tetIslandIds[static_cast<std::size_t>(candidateTet)] = islandId;
                    pendingTets.push(candidateTet);
                }
            }
        }
    }
}
