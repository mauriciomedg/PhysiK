#pragma once

#include <vector>

#include "PhysiK/API/PhysiKAPI.h"
#include "PhysiK/Core/Physics/PhysicsModel.h"
#include "PhysiK/PhysicsData/Node.h"
#include "PhysiK/PhysicsData/Tet.h"

namespace PhysiK
{
    class PHYSIK_API FEMModel : public PhysicsModel
    {
    public:
        void UpdateSystem(World& world, SolverData& solverData, float dt) override;

        static void InitializeTetRestData(Tet& tet, const std::vector<Node>& nodes);
        static void AccumulateElasticForces(
            const std::vector<Tet>& tets,
            const std::vector<Node>& nodes,
            SolverData& solverData);
    };
}
