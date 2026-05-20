#pragma once

#include <string>

#include "PhysiK/API/Handles.h"
#include "PhysiK/API/PhysiKAPI.h"
#include "PhysiK/Components/Component.h"

namespace PhysiK
{
    class PHYSIK_API VisualMeshComponent : public Component
    {
    public:
        VisualMeshComponent();
        VisualMeshComponent(ComponentHandle hostTetMeshHandle, std::string debugEntityName);

        ComponentHandle hostTetMeshHandle;
        std::string debugEntityName;
        bool topologyDirty = false;

        void OnPhysicsEvent(const PhysicsEvent& event) override;
        void Execute(World& world) override;
    };
}
