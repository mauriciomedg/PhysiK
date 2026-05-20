#include "PhysiK/Components/VisualMeshComponent.h"

#include <cstddef>
#include <utility>

#include "PhysiK/Core/World/World.h"

namespace PhysiK
{
    VisualMeshComponent::VisualMeshComponent()
    {
        listenedEvents.push_back(PhysicsEventType::TetMeshTopologyChanged);
    }

    VisualMeshComponent::VisualMeshComponent(
        ComponentHandle hostTetMeshHandle,
        std::string debugEntityName)
        : VisualMeshComponent()
    {
        this->hostTetMeshHandle = hostTetMeshHandle;
        this->debugEntityName = std::move(debugEntityName);
    }

    void VisualMeshComponent::SetVisualMesh(
        const Vec3* vertices,
        int vertexCount,
        const int* triangleIndices,
        int triangleIndexCount)
    {
        restVisualVertices.clear();
        deformedVisualVertices.clear();
        this->triangleIndices.clear();
        embeddedVertices.clear();
        triangleValid.clear();

        if (vertices != nullptr && vertexCount > 0)
        {
            restVisualVertices.assign(vertices, vertices + vertexCount);
            deformedVisualVertices = restVisualVertices;
            embeddedVertices.resize(static_cast<std::size_t>(vertexCount));
        }

        if (triangleIndices != nullptr && triangleIndexCount > 0)
        {
            this->triangleIndices.assign(
                triangleIndices,
                triangleIndices + triangleIndexCount);
            triangleValid.assign(
                static_cast<std::size_t>(triangleIndexCount / 3),
                true);
        }
    }

    const std::vector<Vec3>& VisualMeshComponent::GetDeformedVertices() const
    {
        return deformedVisualVertices;
    }

    const std::vector<int>& VisualMeshComponent::GetTriangleIndices() const
    {
        return triangleIndices;
    }

    void VisualMeshComponent::OnPhysicsEvent(const PhysicsEvent& event)
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

    void VisualMeshComponent::Execute(World& world)
    {
        (void)world;
        if (!topologyDirty)
        {
            return;
        }

        // TODO: rebuild visual mesh embedding/update data after topology changes.
        topologyDirty = false;
    }
}
