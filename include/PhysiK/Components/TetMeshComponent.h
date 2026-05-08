#pragma once

#include <memory>
#include <vector>

#include "PhysiK/Components/Component.h"
#include "PhysiK/Core/Physics/FEM/FEMModel.h"
#include "PhysiK/Math/Vec3.h"
#include "PhysiK/PhysicsData/Material.h"
#include "PhysiK/PhysicsData/Tet.h"

namespace PhysiK
{
    class TetMeshComponent : public Component
    {
    public:
        static std::unique_ptr<TetMeshComponent> CreateFromGlobalNodes(
            World& world,
            const int* globalNodeIndices,
            int nodeCount,
            const int* tetGlobalNodeIndices,
            int tetCount,
            const Material& material);

        static std::unique_ptr<TetMeshComponent> CreateFromPositions(
            World& world,
            const Vec3* positions,
            const float* inverseMasses,
            int nodeCount,
            const int* tetLocalNodeIndices,
            int tetCount,
            const Material& material);

        std::vector<int> nodeIndices;
        std::vector<Tet> tets;

        Material material;
        FEMModel femModel;

        void UpdateSystem(
            World& world,
            SolverData& solverData,
            float dt) override;
    };
}
