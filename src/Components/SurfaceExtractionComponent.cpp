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

        struct TetFaceWithOpposite
        {
            TetFace face;
            int oppositeNode = -1;
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

        TetFaceWithOpposite GetTetFaceWithOpposite(
            const TetMeshComponent& tetMesh,
            int tetIndex,
            int faceIndex)
        {
            const int n0 = tetMesh.GetTetNodeIndex(tetIndex, 0);
            const int n1 = tetMesh.GetTetNodeIndex(tetIndex, 1);
            const int n2 = tetMesh.GetTetNodeIndex(tetIndex, 2);
            const int n3 = tetMesh.GetTetNodeIndex(tetIndex, 3);

            switch (faceIndex)
            {
            case 0:
                return TetFaceWithOpposite{TetFace{n0, n1, n2}, n3};
            case 1:
                return TetFaceWithOpposite{TetFace{n0, n3, n1}, n2};
            case 2:
                return TetFaceWithOpposite{TetFace{n0, n2, n3}, n1};
            case 3:
            default:
                return TetFaceWithOpposite{TetFace{n1, n3, n2}, n0};
            }
        }

        TetFace OrientFaceOutward(
            const TetMeshComponent& tetMesh,
            const TetFace& inputFace,
            int oppositeNode)
        {
            TetFace face = inputFace;
            if (face.node0 < 0 || face.node0 >= tetMesh.GetNodeCount() ||
                face.node1 < 0 || face.node1 >= tetMesh.GetNodeCount() ||
                face.node2 < 0 || face.node2 >= tetMesh.GetNodeCount() ||
                oppositeNode < 0 || oppositeNode >= tetMesh.GetNodeCount())
            {
                return face;
            }

            const Vec3& p0 = tetMesh.GetLocalCurrentPosition(face.node0);
            const Vec3& p1 = tetMesh.GetLocalCurrentPosition(face.node1);
            const Vec3& p2 = tetMesh.GetLocalCurrentPosition(face.node2);
            const Vec3& po = tetMesh.GetLocalCurrentPosition(oppositeNode);

            const Vec3 normal = Cross(p1 - p0, p2 - p0);
            if (normal.LengthSquared() <= 0.000000000001f)
            {
                return face;
            }

            const Vec3 toOpposite = po - p0;
            if (Dot(normal, toOpposite) > 0.0f)
            {
                std::swap(face.node1, face.node2);
            }

            return face;
        }
    }

    ComponentExecutionPriority
    SurfaceExtractionComponent::GetExecutionPriority() const
    {
        return ComponentExecutionPriority::
            SurfaceExtractionComponent;
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
                const TetFaceWithOpposite faceWithOpposite =
                    GetTetFaceWithOpposite(*hostTetMesh, tetIndex, faceIndex);
                const TetFace orientedFace =
                    OrientFaceOutward(
                        *hostTetMesh,
                        faceWithOpposite.face,
                        faceWithOpposite.oppositeNode);
                const FaceKey key = MakeFaceKey(orientedFace);
                FaceEntry& entry = faceEntries[key];
                if (entry.count == 0)
                {
                    entry.orientedFace = orientedFace;
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
