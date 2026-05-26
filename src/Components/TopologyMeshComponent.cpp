#include "PhysiK/Components/TopologyMeshComponent.h"

#include <cstddef>
#include <queue>

#include "PhysiK/Components/TetMeshComponent.h"
#include "PhysiK/Core/World/World.h"
#include "PhysiK/PhysicsData/Tet.h"

namespace PhysiK
{
    namespace
    {
        int CountSharedNodes(const Tet& a, const Tet& b)
        {
            const int aNodes[4] = {a.node0, a.node1, a.node2, a.node3};
            const int bNodes[4] = {b.node0, b.node1, b.node2, b.node3};

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

        bool AreTetNeighbors(const Tet& a, const Tet& b)
        {
            return CountSharedNodes(a, b) >= 3;
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

        const std::vector<Tet>& tets = hostTetMesh->tets;
        tetIslandIds.assign(tets.size(), -1);

        std::queue<int> pendingTets;
        for (int startTet = 0; startTet < static_cast<int>(tets.size()); ++startTet)
        {
            if (!tets[static_cast<std::size_t>(startTet)].active ||
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
                     candidateTet < static_cast<int>(tets.size());
                     ++candidateTet)
                {
                    if (!tets[static_cast<std::size_t>(candidateTet)].active ||
                        tetIslandIds[static_cast<std::size_t>(candidateTet)] != -1)
                    {
                        continue;
                    }

                    if (!AreTetNeighbors(
                            tets[static_cast<std::size_t>(currentTet)],
                            tets[static_cast<std::size_t>(candidateTet)]))
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
