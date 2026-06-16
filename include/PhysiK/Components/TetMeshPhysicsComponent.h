#pragma once

#include <memory>
#include <utility>
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
        struct CachedTetEntry
        {
            int tetIndex = -1;
            int localNodeIndices[4] = {-1, -1, -1, -1};
            int globalNodeIndices[4] = {-1, -1, -1, -1};
            float nodalMass = 0.0f;
            TetFemCache femCache;
            std::pair<int, int> stiffnessNodePairs[16];
        };

        static std::unique_ptr<TetMeshPhysicsComponent> CreateFromGeneratedTetMesh(
            World& world,
            const GeneratedTetMesh& generatedMesh,
            const Material& material);

        static std::unique_ptr<TetMeshPhysicsComponent> CreateFromGeneratedTetMesh(
            World& world,
            const GeneratedTetMesh& generatedMesh,
            const TetMeshPhysicsComponentDesc& desc);

        ComponentExecutionPriority
        GetExecutionPriority() const override;

        int globalNodeBeginIndex = -1;
        int globalNodeCount = 0;
        std::vector<Vec3> nodeVelocities;
        std::vector<TetFemCache> tetFemCache;
        std::vector<CachedTetEntry> cachedTetEntries;
        std::vector<Tet> cachedActiveTets;
        std::vector<TetFemCache> cachedActiveTetFemCache;
        std::vector<int> cachedLocalToGlobalNodeIndices;
        bool femCacheDirty = true;

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

        void MarkFEMCacheDirty();
        void RebuildFEMCacheIfNeeded(const World& world);
        void RebuildFEMCache(const World& world);
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

#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
        int GetFemCacheRebuildCount() const
        {
            return femCacheRebuildCount;
        }

        int GetFemCacheReuseCount() const
        {
            return femCacheReuseCount;
        }
#endif

        void UpdateSystem(
            World& world,
            SolverData& solverData,
            float dt) override;
        void PostUpdate(World& world, float dt) override;

    private:
        int GetCachedGlobalNodeIndex(int localNodeIndex) const;

#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
        int femCacheRebuildCount = 0;
        int femCacheReuseCount = 0;
#endif
    };
}
