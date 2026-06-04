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

    struct TetMeshMappedTet
    {
        int sourceTetIndex = -1;
        bool valid = false;
    };

    class PHYSIK_API TetMeshMapperComponent : public Component
    {
    public:
        TetMeshMapperComponent();
        TetMeshMapperComponent(
            ComponentHandle sourceTetMeshHandle,
            ComponentHandle destinationTetMeshHandle);

        ComponentExecutionPriority
        GetExecutionPriority() const override;

        ComponentHandle sourceTetMeshHandle;
        ComponentHandle destinationTetMeshHandle;
        std::vector<TetMeshMappedVertex> embeddedDestinationVertices;
        std::vector<TetMeshMappedTet> embeddedDestinationTets;

        void MarkMappingDirty();
        void MarkActiveStateDirty();
        bool IsMappingDirty() const;
        void OnPhysicsEvent(const PhysicsEvent& event) override;
        void PostUpdate(World& world, float dt) override;

    private:
        bool BuildTetMeshMapping(World& world);
        bool RefreshDestinationActiveStates(World& world);
        void UpdateDestinationNodes(World& world);

        bool mappingDirty = true;
        bool activeStateDirty = true;
    };
}
