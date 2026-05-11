#pragma once

#include <cstdint>
#include <vector>

#include "PhysiK/API/PhysiKAPI.h"
#include "PhysiK/Core/Physics/PhysicsModel.h"
#include "PhysiK/PhysicsData/Node.h"
#include "PhysiK/PhysicsData/Tet.h"

namespace PhysiK
{
    enum class FemModel : std::uint32_t
    {
        Linear = 0,
        Corotational = 1,
        NeoHookean = 2
    };

    class SolverData;
    class TetMeshComponent;
    class World;

    class PHYSIK_API FEMModel : public PhysicsModel
    {
    public:
        void UpdateSystem(
            World& world,
            TetMeshComponent& owner,
            SolverData& solverData,
            float dt);

        static void InitializeTetRestData(Tet& tet, const std::vector<Node>& nodes);
        static void AccumulateElasticForces(
            const std::vector<Tet>& tets,
            const std::vector<Node>& nodes,
            SolverData& solverData);
        static void AccumulateCorotationalElasticForces(
            const std::vector<Tet>& tets,
            const std::vector<Node>& nodes,
            SolverData& solverData);
        static bool AccumulateForces(
            FemModel femModel,
            const std::vector<Tet>& tets,
            const std::vector<Node>& nodes,
            SolverData& solverData);
        static bool IsFemModelImplemented(FemModel femModel);
        static const char* GetNotImplementedMessage(FemModel femModel);
    };
}
