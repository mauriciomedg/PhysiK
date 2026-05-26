#pragma once

#include <vector>

#include "PhysiK/API/Handles.h"
#include "PhysiK/API/PhysiKAPI.h"
#include "PhysiK/Components/Component.h"
#include "PhysiK/Math/Vec3.h"
#include "PhysiK/Math/Vec4.h"

namespace PhysiK
{
    struct VisualTetEmbeddedVertex
    {
        int mechanicalTetIndex = -1;
        Vec4 barycentric;
        bool valid = false;
    };

    class PHYSIK_API VisualTetMeshComponent : public Component
    {
    public:
        VisualTetMeshComponent() = default;
        explicit VisualTetMeshComponent(ComponentHandle hostMechanicalTetMeshHandle);

        ComponentHandle hostMechanicalTetMeshHandle;
        std::vector<Vec3> restVisualVertices;
        std::vector<Vec3> deformedVisualVertices;
        std::vector<int> visualTetIndices;
        std::vector<VisualTetEmbeddedVertex> embeddedVertices;
        std::vector<bool> activeVisualTets;

        void SetVisualTetMeshData(
            const Vec3* vertices,
            int vertexCount,
            const int* tetIndices,
            int tetIndexCount);
        void BuildEmbedding(const World& world);
        void UpdateDeformedVertices(const World& world);
        const std::vector<Vec3>& GetDeformedVertices() const;
        const std::vector<int>& GetVisualTetIndices() const;
        void PostUpdate(World& world, float dt) override;
    };
}
