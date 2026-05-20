#include "PhysiK/Components/VisualMeshComponent.h"

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
