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
        FemModel femModel = FemModel::Corotational;
    };

    class PHYSIK_API TetMeshPhysicsComponent : public TetMeshComponent
    {
    public:
        static std::unique_ptr<TetMeshPhysicsComponent> CreateFromGeneratedTetMesh(
            World& world,
            const GeneratedTetMesh& generatedMesh,
            const Material& material);

        static std::unique_ptr<TetMeshPhysicsComponent> CreateFromGeneratedTetMesh(
            World& world,
            const GeneratedTetMesh& generatedMesh,
            const TetMeshPhysicsComponentDesc& desc);

        int globalNodeBeginIndex = -1;
        int globalNodeCount = 0;
        std::vector<Vec3> nodeVelocities;
        std::vector<TetFemCache> tetFemCache;

        Material material;
        FemModel selectedFemModel = FemModel::Corotational;
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

        void RebuildFemRestData();
        void RebuildTetFemCache();
        void EnsureFemSparsePattern(int worldNodeCount);
        void SyncCurrentPositionsFromWorld(const World& world);

        int GetGlobalNodeBeginIndex() const;
        int GetGlobalNodeCount() const;
        int GetGlobalNodeIndex(int localNodeIndex) const override;
        void SetLocalCurrentPosition(int localNodeIndex, const Vec3& position) override;
        bool SetTetActive(int tetIndex, bool active) override;
        bool DeactivateTet(int tetIndex) override;

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
