#pragma once

#include <vector>

#include "PhysiK/API/Handles.h"
#include "PhysiK/API/PhysiKAPI.h"
#include "PhysiK/Components/Component.h"

namespace PhysiK
{
    class PHYSIK_API TopologyMeshComponent : public Component
    {
    public:
        TopologyMeshComponent();
        explicit TopologyMeshComponent(ComponentHandle hostTetMeshHandle);

        ComponentExecutionPriority
        GetExecutionPriority() const override;

        ComponentHandle hostTetMeshHandle;
        bool topologyDirty = true;
        std::vector<int> tetIslandIds;
        int islandCount = 0;

        int GetTetIslandId(int tetIndex) const;
        int GetIslandCount() const;

        void OnPhysicsEvent(const PhysicsEvent& event) override;
        void PostUpdate(World& world, float dt) override;

    private:
        void RebuildTopology(World& world);
    };
}
