#pragma once

#include <memory>
#include <vector>

#include "PhysiK/Components/TetMeshComponent.h"
#include "PhysiK/Core/Physics/FEM/FEMModel.h"
#include "PhysiK/Math/SparseBlockMatrix.h"
#include "PhysiK/PhysicsData/Material.h"

namespace PhysiK
{
    struct TetMeshPhysicsComponentDesc
    {
        Material material;
        FemModel femModel = FemModel::Linear;
    };

    class PHYSIK_API TetMeshPhysicsComponent : public TetMeshComponent
    {
    public:
        static std::unique_ptr<TetMeshPhysicsComponent> CreateFromGlobalNodes(
            World& world,
            const int* globalNodeIndices,
            int nodeCount,
            const int* tetGlobalNodeIndices,
            int tetCount,
            const Material& material);

        static std::unique_ptr<TetMeshPhysicsComponent> CreateFromGlobalNodes(
            World& world,
            const int* globalNodeIndices,
            int nodeCount,
            const int* tetGlobalNodeIndices,
            int tetCount,
            const TetMeshPhysicsComponentDesc& desc);

        static std::unique_ptr<TetMeshPhysicsComponent> CreateFromPositions(
            World& world,
            const Vec3* positions,
            const int* fixedNodeFlags,
            int nodeCount,
            const int* tetLocalNodeIndices,
            int tetCount,
            const Material& material);

        static std::unique_ptr<TetMeshPhysicsComponent> CreateFromPositions(
            World& world,
            const Vec3* positions,
            const int* fixedNodeFlags,
            int nodeCount,
            const int* tetLocalNodeIndices,
            int tetCount,
            const TetMeshPhysicsComponentDesc& desc);

        std::vector<int> worldNodeIndices;
        std::vector<Tet> worldTets;
        std::vector<TetFemCache> tetFemCache;

        Material material;
        FemModel selectedFemModel = FemModel::Linear;
        FEMModel femModel;
        SparseBlockMatrix femSparseMatrix;
        bool femSparsePatternDirty = true;
        bool physicsEnabled = true;

        void SetMaterial(const Material& value);
        void SetFemModel(FemModel value)
        {
            selectedFemModel = value;
        }

        FemModel GetFemModel() const
        {
            return selectedFemModel;
        }

        void MarkFemSparsePatternDirty()
        {
            femSparsePatternDirty = true;
        }

        void RebuildWorldTets(const World& world);
        void RebuildTetFemCache();
        void EnsureFemSparsePattern(int worldNodeCount);
        void SyncWorldTetActiveStates();
        void SyncCurrentPositionsFromWorld(const World& world);

        int GetWorldNodeIndex(int localNodeIndex) const override;
        void SetLocalCurrentPosition(int localNodeIndex, const Vec3& position) override;
        void SetTetActive(int tetIndex, bool active) override;
        void DeactivateTet(int tetIndex) override;

        const SparseBlockMatrix& GetFemSparseMatrix() const
        {
            return femSparseMatrix;
        }

        void UpdateSystem(
            World& world,
            SolverData& solverData,
            float dt) override;
        void PostUpdate(World& world, float dt) override;
    };
}
