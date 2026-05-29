#include "PhysiK/Components/SurfaceExtractionComponent.h"

#include <cstddef>
#include <functional>
#include <unordered_map>
#include <utility>

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

        struct FaceEntry
        {
            TetFace orientedFace;
            int count = 0;
        };

        struct FaceKey
        {
            int node0 = -1;
            int node1 = -1;
            int node2 = -1;

            bool operator==(const FaceKey& other) const
            {
                return node0 == other.node0 &&
                    node1 == other.node1 &&
                    node2 == other.node2;
            }
        };

        struct FaceKeyHash
        {
            std::size_t operator()(const FaceKey& key) const
            {
                const std::size_t h0 = std::hash<int>{}(key.node0);
                const std::size_t h1 = std::hash<int>{}(key.node1);
                const std::size_t h2 = std::hash<int>{}(key.node2);
                return h0 ^ (h1 << 1u) ^ (h2 << 2u);
            }
        };

        FaceKey MakeFaceKey(const TetFace& face)
        {
            FaceKey key{face.node0, face.node1, face.node2};
            if (key.node1 < key.node0)
            {
                std::swap(key.node0, key.node1);
            }
            if (key.node2 < key.node1)
            {
                std::swap(key.node1, key.node2);
            }
            if (key.node1 < key.node0)
            {
                std::swap(key.node0, key.node1);
            }
            return key;
        }

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

        std::unordered_map<FaceKey, FaceEntry, FaceKeyHash> faceEntries;
        for (int tetIndex = 0; tetIndex < hostTetMesh->GetTetCount(); ++tetIndex)
        {
            if (!hostTetMesh->IsTetActive(tetIndex))
            {
                continue;
            }

            for (int faceIndex = 0; faceIndex < 4; ++faceIndex)
            {
                const TetFace face = GetTetFace(*hostTetMesh, tetIndex, faceIndex);
                const FaceKey key = MakeFaceKey(face);
                FaceEntry& entry = faceEntries[key];
                if (entry.count == 0)
                {
                    entry.orientedFace = face;
                }

                ++entry.count;
            }
        }

        for (const auto& item : faceEntries)
        {
            const FaceEntry& entry = item.second;
            if (entry.count != 1)
            {
                continue;
            }

            surfaceTriangleIndices.push_back(entry.orientedFace.node0);
            surfaceTriangleIndices.push_back(entry.orientedFace.node1);
            surfaceTriangleIndices.push_back(entry.orientedFace.node2);
        }
    }

    const std::vector<int>& SurfaceExtractionComponent::GetSurfaceTriangleIndices() const
    {
        return surfaceTriangleIndices;
    }

    ComponentHandle SurfaceExtractionComponent::GetHostTetMeshHandle() const
    {
        return hostTetMeshHandle;
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
