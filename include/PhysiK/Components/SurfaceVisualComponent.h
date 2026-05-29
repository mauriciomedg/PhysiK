#pragma once

#include <vector>

#include "PhysiK/API/Handles.h"
#include "PhysiK/API/PhysiKAPI.h"
#include "PhysiK/Components/Component.h"
#include "PhysiK/Math/Vec3.h"

namespace PhysiK
{
    class World;

    class PHYSIK_API SurfaceVisualComponent final : public Component
    {
    public:
        SurfaceVisualComponent();
        explicit SurfaceVisualComponent(ComponentHandle surfaceExtractionHandle);

        void RebuildVisualSurface(const World& world);

        const std::vector<Vec3>& GetVisualVertices() const;
        const std::vector<int>& GetVisualTriangleIndices() const;
        const std::vector<Vec3>& GetVisualNormals() const;

        ComponentHandle GetSurfaceExtractionHandle() const;

        void OnPhysicsEvent(const PhysicsEvent& event) override;
        void PostUpdate(World& world, float dt) override;

    private:
        ComponentHandle surfaceExtractionHandle;

        std::vector<Vec3> visualVertices;
        std::vector<int> visualTriangleIndices;
        std::vector<Vec3> visualNormals;

        bool visualDirty = true;
    };
}
