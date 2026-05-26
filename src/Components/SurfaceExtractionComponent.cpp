#include "PhysiK/Components/SurfaceExtractionComponent.h"

#include <cstddef>

#include "PhysiK/Components/TetMeshComponent.h"
#include "PhysiK/Core/World/World.h"
#include "PhysiK/PhysicsData/Tet.h"

namespace PhysiK
{
    namespace
    {
        struct TetFace
        {
            int node0 = -1;
            int node1 = -1;
            int node2 = -1;
        };

        TetFace GetTetFace(const Tet& tet, int faceIndex)
        {
            switch (faceIndex)
            {
            case 0:
                return TetFace{tet.node0, tet.node1, tet.node2};
            case 1:
                return TetFace{tet.node0, tet.node3, tet.node1};
            case 2:
                return TetFace{tet.node0, tet.node2, tet.node3};
            case 3:
            default:
                return TetFace{tet.node1, tet.node3, tet.node2};
            }
        }

        bool HasSameUndirectedNodes(const TetFace& a, const TetFace& b)
        {
            const int aNodes[3] = {a.node0, a.node1, a.node2};
            const int bNodes[3] = {b.node0, b.node1, b.node2};

            for (int aNode : aNodes)
            {
                bool found = false;
                for (int bNode : bNodes)
                {
                    if (aNode == bNode)
                    {
                        found = true;
                        break;
                    }
                }

                if (!found)
                {
                    return false;
                }
            }

            return true;
        }

        bool HasActiveNeighborSharingFace(
            const std::vector<Tet>& tets,
            int sourceTetIndex,
            const TetFace& face)
        {
            for (int tetIndex = 0; tetIndex < static_cast<int>(tets.size()); ++tetIndex)
            {
                if (tetIndex == sourceTetIndex ||
                    !tets[static_cast<std::size_t>(tetIndex)].active)
                {
                    continue;
                }

                const Tet& candidateTet = tets[static_cast<std::size_t>(tetIndex)];
                for (int faceIndex = 0; faceIndex < 4; ++faceIndex)
                {
                    if (HasSameUndirectedNodes(face, GetTetFace(candidateTet, faceIndex)))
                    {
                        return true;
                    }
                }
            }

            return false;
        }
    }

    SurfaceExtractionComponent::SurfaceExtractionComponent()
    {
        listenedEvents.push_back(PhysicsEventType::TetMeshTopologyChanged);
    }

    SurfaceExtractionComponent::SurfaceExtractionComponent(ComponentHandle hostTetMeshHandle)
        : SurfaceExtractionComponent()
    {
        this->hostTetMeshHandle = hostTetMeshHandle;
    }

    void SurfaceExtractionComponent::RebuildSurface(const World& world)
    {
        surfaceTriangleIndices.clear();

        const Component* hostComponent = world.GetComponent(hostTetMeshHandle);
        const TetMeshComponent* hostTetMesh =
            dynamic_cast<const TetMeshComponent*>(hostComponent);
        if (hostTetMesh == nullptr)
        {
            return;
        }

        const std::vector<Tet>& tets = hostTetMesh->tets;
        surfaceTriangleIndices.reserve(tets.size() * 12u);

        for (int tetIndex = 0; tetIndex < static_cast<int>(tets.size()); ++tetIndex)
        {
            const Tet& tet = tets[static_cast<std::size_t>(tetIndex)];
            if (!tet.active)
            {
                continue;
            }

            for (int faceIndex = 0; faceIndex < 4; ++faceIndex)
            {
                const TetFace face = GetTetFace(tet, faceIndex);
                if (HasActiveNeighborSharingFace(tets, tetIndex, face))
                {
                    continue;
                }

                surfaceTriangleIndices.push_back(face.node0);
                surfaceTriangleIndices.push_back(face.node1);
                surfaceTriangleIndices.push_back(face.node2);
            }
        }
    }

    const std::vector<int>& SurfaceExtractionComponent::GetSurfaceTriangleIndices() const
    {
        return surfaceTriangleIndices;
    }

    void SurfaceExtractionComponent::OnPhysicsEvent(const PhysicsEvent& event)
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

        surfaceDirty = true;
    }

    void SurfaceExtractionComponent::PostUpdate(World& world, float dt)
    {
        (void)dt;

        if (!surfaceDirty)
        {
            return;
        }

        RebuildSurface(world);
        surfaceDirty = false;
    }
}
