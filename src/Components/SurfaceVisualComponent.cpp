#include "PhysiK/Components/SurfaceVisualComponent.h"

#include <cstddef>

#include "PhysiK/Components/SurfaceExtractionComponent.h"
#include "PhysiK/Components/TetMeshComponent.h"
#include "PhysiK/Core/World/World.h"

namespace PhysiK
{
    namespace
    {
        Vec3 NormalizeOrZero(const Vec3& vector)
        {
            const float length = vector.Length();
            if (length <= 0.0f)
            {
                return Vec3{};
            }

            return vector / length;
        }
    }

    SurfaceVisualComponent::SurfaceVisualComponent()
    {
        listenedEvents.push_back(PhysicsEventType::TetMeshTopologyChanged);
    }

    SurfaceVisualComponent::SurfaceVisualComponent(
        ComponentHandle surfaceExtractionHandle)
        : SurfaceVisualComponent()
    {
        this->surfaceExtractionHandle = surfaceExtractionHandle;
    }

    void SurfaceVisualComponent::RebuildVisualSurface(const World& world)
    {
        visualVertices.clear();
        visualTriangleIndices.clear();
        visualNormals.clear();

        const Component* surfaceComponent =
            world.GetComponent(surfaceExtractionHandle);
        const SurfaceExtractionComponent* surfaceExtraction =
            dynamic_cast<const SurfaceExtractionComponent*>(surfaceComponent);
        if (surfaceExtraction == nullptr)
        {
            return;
        }

        const Component* hostComponent =
            world.GetComponent(surfaceExtraction->GetHostTetMeshHandle());
        const TetMeshComponent* hostTetMesh =
            dynamic_cast<const TetMeshComponent*>(hostComponent);
        if (hostTetMesh == nullptr)
        {
            return;
        }

        const std::vector<int>& rawTriangles =
            surfaceExtraction->GetSurfaceTriangleIndices();
        const std::size_t triangleCount = rawTriangles.size() / 3u;
        visualVertices.reserve(triangleCount * 3u);
        visualTriangleIndices.reserve(triangleCount * 3u);
        visualNormals.reserve(triangleCount * 3u);

        for (std::size_t triangleIndex = 0;
             triangleIndex < triangleCount;
             ++triangleIndex)
        {
            const std::size_t rawIndex = triangleIndex * 3u;
            const int node0 = rawTriangles[rawIndex + 0u];
            const int node1 = rawTriangles[rawIndex + 1u];
            const int node2 = rawTriangles[rawIndex + 2u];

            if (node0 < 0 || node0 >= hostTetMesh->GetNodeCount() ||
                node1 < 0 || node1 >= hostTetMesh->GetNodeCount() ||
                node2 < 0 || node2 >= hostTetMesh->GetNodeCount())
            {
                continue;
            }

            const Vec3& p0 = hostTetMesh->GetLocalCurrentPosition(node0);
            const Vec3& p1 = hostTetMesh->GetLocalCurrentPosition(node1);
            const Vec3& p2 = hostTetMesh->GetLocalCurrentPosition(node2);
            const Vec3 normal = NormalizeOrZero(Cross(p1 - p0, p2 - p0));

            const int baseIndex = static_cast<int>(visualVertices.size());
            visualVertices.push_back(p0);
            visualVertices.push_back(p1);
            visualVertices.push_back(p2);

            visualNormals.push_back(normal);
            visualNormals.push_back(normal);
            visualNormals.push_back(normal);

            visualTriangleIndices.push_back(baseIndex + 0);
            visualTriangleIndices.push_back(baseIndex + 1);
            visualTriangleIndices.push_back(baseIndex + 2);
        }
    }

    const std::vector<Vec3>& SurfaceVisualComponent::GetVisualVertices() const
    {
        return visualVertices;
    }

    const std::vector<int>& SurfaceVisualComponent::GetVisualTriangleIndices() const
    {
        return visualTriangleIndices;
    }

    const std::vector<Vec3>& SurfaceVisualComponent::GetVisualNormals() const
    {
        return visualNormals;
    }

    ComponentHandle SurfaceVisualComponent::GetSurfaceExtractionHandle() const
    {
        return surfaceExtractionHandle;
    }

    void SurfaceVisualComponent::OnPhysicsEvent(const PhysicsEvent& event)
    {
        if (event.type != PhysicsEventType::TetMeshTopologyChanged ||
            event.world == nullptr)
        {
            return;
        }

        const Component* surfaceComponent =
            event.world->GetComponent(surfaceExtractionHandle);
        const SurfaceExtractionComponent* surfaceExtraction =
            dynamic_cast<const SurfaceExtractionComponent*>(surfaceComponent);
        if (surfaceExtraction == nullptr)
        {
            return;
        }

        if (event.world->GetComponent(surfaceExtraction->GetHostTetMeshHandle()) !=
            event.sender)
        {
            return;
        }

        visualDirty = true;
    }

    void SurfaceVisualComponent::PostUpdate(World& world, float dt)
    {
        (void)dt;

        if (!visualDirty)
        {
            return;
        }

        RebuildVisualSurface(world);
        visualDirty = false;
    }
}
