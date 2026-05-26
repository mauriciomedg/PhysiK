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

        TetFace GetTetFace(const TetMeshComponent& tetMesh, int tetIndex, int faceIndex)
        {
            switch (faceIndex)
            {
            case 0:
                return TetFace{
                    tetMesh.GetTetNodeIndex(tetIndex, 0),
                    tetMesh.GetTetNodeIndex(tetIndex, 1),
                    tetMesh.GetTetNodeIndex(tetIndex, 2)};
            case 1:
                return TetFace{
                    tetMesh.GetTetNodeIndex(tetIndex, 0),
                    tetMesh.GetTetNodeIndex(tetIndex, 3),
                    tetMesh.GetTetNodeIndex(tetIndex, 1)};
            case 2:
                return TetFace{
                    tetMesh.GetTetNodeIndex(tetIndex, 0),
                    tetMesh.GetTetNodeIndex(tetIndex, 2),
                    tetMesh.GetTetNodeIndex(tetIndex, 3)};
            case 3:
            default:
                return TetFace{
                    tetMesh.GetTetNodeIndex(tetIndex, 1),
                    tetMesh.GetTetNodeIndex(tetIndex, 3),
                    tetMesh.GetTetNodeIndex(tetIndex, 2)};
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
            const TetMeshComponent& tetMesh,
            int sourceTetIndex,
            const TetFace& face)
        {
            for (int tetIndex = 0; tetIndex < tetMesh.GetTetCount(); ++tetIndex)
            {
                if (tetIndex == sourceTetIndex ||
                    !tetMesh.IsTetActive(tetIndex))
                {
                    continue;
                }

                for (int faceIndex = 0; faceIndex < 4; ++faceIndex)
                {
                    if (HasSameUndirectedNodes(
                            face,
                            GetTetFace(tetMesh, tetIndex, faceIndex)))
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

        surfaceTriangleIndices.reserve(
            static_cast<std::size_t>(hostTetMesh->GetTetCount()) * 12u);

        for (int tetIndex = 0; tetIndex < hostTetMesh->GetTetCount(); ++tetIndex)
        {
            if (!hostTetMesh->IsTetActive(tetIndex))
            {
                continue;
            }

            for (int faceIndex = 0; faceIndex < 4; ++faceIndex)
            {
                const TetFace face = GetTetFace(*hostTetMesh, tetIndex, faceIndex);
                if (HasActiveNeighborSharingFace(*hostTetMesh, tetIndex, face))
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
