#pragma once

#include <vector>

#include "PhysiK/API/Handles.h"
#include "PhysiK/API/PhysiKAPI.h"
#include "PhysiK/Components/Component.h"
#include "PhysiK/Math/Vec4.h"

namespace PhysiK
{
    struct TetMeshMappedVertex
    {
        int sourceTetIndex = -1;
        Vec4 barycentric;
        bool valid = false;
    };

    class PHYSIK_API TetMeshMapperComponent : public Component
    {
    public:
        TetMeshMapperComponent() = default;
        TetMeshMapperComponent(
            ComponentHandle sourceTetMeshHandle,
            ComponentHandle destinationTetMeshHandle);

        ComponentHandle sourceTetMeshHandle;
        ComponentHandle destinationTetMeshHandle;
        std::vector<TetMeshMappedVertex> embeddedDestinationVertices;
        bool mappingDirty = true;

        bool BuildTetMeshMapping(World& world);
        void UpdateDestinationNodes(World& world);
        void MarkMappingDirty();
        bool IsMappingDirty() const;
        void PostUpdate(World& world, float dt) override;
    };
}
