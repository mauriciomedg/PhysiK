#pragma once

#include <array>
#include <cstdint>
#include <vector>

#include "PhysiK/API/PhysiKAPI.h"
#include "PhysiK/Core/Physics/PhysicsModel.h"
#include "PhysiK/PhysicsData/Material.h"
#include "PhysiK/PhysicsData/Node.h"
#include "PhysiK/PhysicsData/Tet.h"

namespace PhysiK
{
    using Matrix6 = std::array<std::array<float, 6>, 6>;
    using Matrix6x12 = std::array<std::array<float, 12>, 6>;
    using Matrix12 = std::array<std::array<float, 12>, 12>;

    struct TetFemCache
    {
        Matrix6x12 B{};
        Matrix6 D{};
        Matrix12 Ke{};
    };

    enum class FemModel : std::uint32_t
    {
        Linear = 0,
        Corotational = 1,
        NeoHookean = 2
    };

    class SolverData;
    class TetMeshPhysicsComponent;
    class World;

    class PHYSIK_API FEMModel : public PhysicsModel
    {
    public:
        void UpdateSystem(
            World& world,
            TetMeshPhysicsComponent& owner,
            SolverData& solverData,
            float dt);

        static void InitializeTetRestData(Tet& tet, const std::vector<Node>& nodes);
        static TetFemCache BuildTetFemCache(const Tet& tet);
        static void AccumulateElasticForces(
            const std::vector<Tet>& tets,
            const std::vector<Node>& nodes,
            SolverData& solverData);
        static void AccumulateElasticForces(
            const std::vector<Tet>& tets,
            const std::vector<TetFemCache>& tetFemCache,
            const std::vector<Node>& nodes,
            SolverData& solverData);
        static void AccumulateCorotationalElasticForces(
            const std::vector<Tet>& tets,
            const std::vector<Node>& nodes,
            SolverData& solverData);
        static void AccumulateCorotationalElasticForces(
            const std::vector<Tet>& tets,
            const std::vector<TetFemCache>& tetFemCache,
            const std::vector<Node>& nodes,
            SolverData& solverData);
        static bool AccumulateForces(
            FemModel femModel,
            const std::vector<Tet>& tets,
            const std::vector<Node>& nodes,
            SolverData& solverData);
        static bool AccumulateForces(
            FemModel femModel,
            const std::vector<Tet>& tets,
            const std::vector<TetFemCache>& tetFemCache,
            const std::vector<Node>& nodes,
            SolverData& solverData);

        static void AddLumpedMassToSolverData(
            const World& world,
            SolverData& solverData,
            int nodeIndex,
            float mass);

        static void AssembleLumpedMass(
            const TetMeshPhysicsComponent& component,
            const World& world,
            SolverData& solverData);
        static void AssembleLumpedMass(
            const Material& material,
            const std::vector<Tet>& tets,
            const World& world,
            SolverData& solverData);
       
        static std::vector<std::pair<int, int>> BuildSparsePatternFromTetConnectivity(
            const std::vector<Tet>& tets);

        static bool IsFemModelImplemented(FemModel femModel);
        static const char* GetNotImplementedMessage(FemModel femModel);
    };
}
