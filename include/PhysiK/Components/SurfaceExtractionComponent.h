#pragma once

#include <vector>

#include "PhysiK/API/Handles.h"
#include "PhysiK/API/PhysiKAPI.h"
#include "PhysiK/Components/Component.h"

namespace PhysiK
{
    class PHYSIK_API SurfaceExtractionComponent : public Component
    {
    public:
        SurfaceExtractionComponent();
        explicit SurfaceExtractionComponent(ComponentHandle hostTetMeshHandle);

        ComponentExecutionPriority
        GetExecutionPriority() const override;

        ComponentHandle hostTetMeshHandle;
        std::vector<int> surfaceTriangleIndices;
        bool surfaceDirty = true;

        void RebuildSurface(const World& world);
        const std::vector<int>& GetSurfaceTriangleIndices() const;
        ComponentHandle GetHostTetMeshHandle() const;

        void OnPhysicsEvent(const PhysicsEvent& event) override;
        void PostUpdate(World& world, float dt) override;
    };
}
