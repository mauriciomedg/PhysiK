#pragma once

#include <memory>
#include <vector>

#include "PhysiK/API/PhysiKAPI.h"
#include "PhysiK/Components/Component.h"
#include "PhysiK/Core/Physics/FEM/FEMModel.h"
#include "PhysiK/Math/SparseBlockMatrix.h"
#include "PhysiK/Math/Vec3.h"
#include "PhysiK/PhysicsData/Material.h"
#include "PhysiK/PhysicsData/Tet.h"

namespace PhysiK
{
    struct TetMeshComponentDesc
    {
        Material material;
        FemModel femModel = FemModel::Linear;
    };

    class PHYSIK_API TetMeshComponent : public Component
    {
    public:
        static std::unique_ptr<TetMeshComponent> CreateFromGlobalNodes(
            World& world,
            const int* globalNodeIndices,
            int nodeCount,
            const int* tetGlobalNodeIndices,
            int tetCount,
            const Material& material);

        static std::unique_ptr<TetMeshComponent> CreateFromGlobalNodes(
            World& world,
            const int* globalNodeIndices,
            int nodeCount,
            const int* tetGlobalNodeIndices,
            int tetCount,
            const TetMeshComponentDesc& desc);

        static std::unique_ptr<TetMeshComponent> CreateFromPositions(
            World& world,
            const Vec3* positions,
            const int* fixedNodeFlags,
            int nodeCount,
            const int* tetLocalNodeIndices,
            int tetCount,
            const Material& material);

        static std::unique_ptr<TetMeshComponent> CreateFromPositions(
            World& world,
            const Vec3* positions,
            const int* fixedNodeFlags,
            int nodeCount,
            const int* tetLocalNodeIndices,
            int tetCount,
            const TetMeshComponentDesc& desc);

        std::vector<int> nodeIndices;
        std::vector<Tet> tets;
        std::vector<TetFemCache> tetFemCache;

        Material material;
        FemModel selectedFemModel = FemModel::Linear;
        FEMModel femModel;
        SparseBlockMatrix femSparseMatrix;
        bool femSparsePatternDirty = true;

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

        void RebuildTetFemCache();
        void EnsureFemSparsePattern(int worldNodeCount);

        bool IsTetActive(int tetIndex) const;
        void SetTetActive(int tetIndex, bool active);
        void DeactivateTet(int tetIndex);
        int GetActiveTetCount() const;

        const SparseBlockMatrix& GetFemSparseMatrix() const
        {
            return femSparseMatrix;
        }

        void UpdateSystem(
            World& world,
            SolverData& solverData,
            float dt) override;
    };
}
