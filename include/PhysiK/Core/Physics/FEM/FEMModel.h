#pragma once

#include <array>
#include <cstdint>
#include <vector>

#include "PhysiK/API/PhysiKAPI.h"
#include "PhysiK/Core/Physics/PhysicsModel.h"
#include "PhysiK/Math/Mat3.h"
#include "PhysiK/Math/Vec3.h"
#include "PhysiK/PhysicsData/Material.h"
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

    struct TetElementContribution
    {
        int localNodeIndices[4] = {-1, -1, -1, -1};
        Vec3 forces[4];
        Mat3 stiffness[4][4];
    };

    struct TetMassContribution
    {
        int localNodeIndices[4] = {-1, -1, -1, -1};
        float nodalMass = 0.0f;
    };

    enum class FemModel : std::uint32_t
    {
        Linear = 0,
        Corotational = 1,
        NeoHookean = 2
    };

    class PHYSIK_API FEMModel : public PhysicsModel
    {
    public:
        static void InitializeTetRestData(
            Tet& tet,
            const std::vector<Vec3>& restPositions);
        static TetFemCache BuildTetFemCache(const Tet& tet);
        static void ComputeElasticForces(
            const std::vector<Tet>& tets,
            const std::vector<Vec3>& positions,
            const std::vector<Vec3>& velocities,
            std::vector<TetElementContribution>& outContributions);
        static void ComputeElasticForces(
            const std::vector<Tet>& tets,
            const std::vector<TetFemCache>& tetFemCache,
            const std::vector<Vec3>& positions,
            const std::vector<Vec3>& velocities,
            std::vector<TetElementContribution>& outContributions);
        static void ComputeCorotationalElasticForces(
            const std::vector<Tet>& tets,
            const std::vector<Vec3>& positions,
            const std::vector<Vec3>& velocities,
            std::vector<TetElementContribution>& outContributions);
        static void ComputeCorotationalElasticForces(
            const std::vector<Tet>& tets,
            const std::vector<TetFemCache>& tetFemCache,
            const std::vector<Vec3>& positions,
            const std::vector<Vec3>& velocities,
            std::vector<TetElementContribution>& outContributions);
        static bool ComputeForces(
            FemModel femModel,
            const std::vector<Tet>& tets,
            const std::vector<Vec3>& positions,
            const std::vector<Vec3>& velocities,
            std::vector<TetElementContribution>& outContributions);
        static bool ComputeForces(
            FemModel femModel,
            const std::vector<Tet>& tets,
            const std::vector<TetFemCache>& tetFemCache,
            const std::vector<Vec3>& positions,
            const std::vector<Vec3>& velocities,
            std::vector<TetElementContribution>& outContributions);

        static void ComputeLumpedMass(
            const Material& material,
            const std::vector<Tet>& tets,
            std::vector<TetMassContribution>& outContributions);
       
        static std::vector<std::pair<int, int>> BuildSparsePatternFromTetConnectivity(
            const std::vector<Tet>& tets);

        static bool IsFemModelImplemented(FemModel femModel);
        static const char* GetNotImplementedMessage(FemModel femModel);
    };
}
