#pragma once

#include <string>
#include <vector>

#include "PhysiK/API/Handles.h"
#include "PhysiK/API/PhysiKAPI.h"
#include "PhysiK/Components/Component.h"
#include "PhysiK/Math/Vec3.h"
#include "PhysiK/Math/Vec4.h"

namespace PhysiK
{
    struct EmbeddedVertex
    {
        int tetIndex = -1;
        Vec4 barycentric;
        bool valid = false;
    };

    class PHYSIK_API VisualMeshComponent : public Component
    {
    public:
        VisualMeshComponent();
        VisualMeshComponent(ComponentHandle hostTetMeshHandle, std::string debugEntityName);

        ComponentHandle hostTetMeshHandle;
        std::string debugEntityName;
        bool topologyDirty = false;
        std::vector<Vec3> restVisualVertices;
        std::vector<Vec3> deformedVisualVertices;
        std::vector<int> triangleIndices;
        std::vector<EmbeddedVertex> embeddedVertices;
        std::vector<bool> triangleValid;

        void SetVisualMesh(
            const Vec3* vertices,
            int vertexCount,
            const int* triangleIndices,
            int triangleIndexCount);
        const std::vector<Vec3>& GetDeformedVertices() const;
        const std::vector<int>& GetTriangleIndices() const;
        void OnPhysicsEvent(const PhysicsEvent& event) override;
        void Execute(World& world) override;
    };
}
