#pragma once

#include <vector>

#include "PhysiK/Components/Component.h"
#include "PhysiK/Core/Physics/FEM/FEMModel.h"
#include "PhysiK/PhysicsData/Material.h"
#include "PhysiK/PhysicsData/Tet.h"

namespace PhysiK
{
    class TetMeshComponent : public Component
    {
    public:
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
