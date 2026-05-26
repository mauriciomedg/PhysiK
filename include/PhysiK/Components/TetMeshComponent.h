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
        std::vector<Vec3> restNodePositions;
        std::vector<Vec3> currentNodePositions;
        std::vector<int> nodeIndices;
        std::vector<Tet> tets;

        bool topologyDirty = false;

        void SetGeometry(
            const Vec3* positions,
            int nodeCount,
            const int* tetLocalNodeIndices,
            int tetCount);

        int GetNodeCount() const;
        int GetTetCount() const;
        int GetTetNodeIndex(int tetIndex, int cornerIndex) const;
        const Vec3& GetLocalRestPosition(int localNodeIndex) const;
        const Vec3& GetLocalCurrentPosition(int localNodeIndex) const;
        virtual int GetWorldNodeIndex(int localNodeIndex) const;
        virtual void SetLocalCurrentPosition(int localNodeIndex, const Vec3& position);

        bool IsTetActive(int tetIndex) const;
        void SetTetActive(int tetIndex, bool active);
        void DeactivateTet(int tetIndex);
        int GetActiveTetCount() const;

        void PostUpdate(World& world, float dt) override;
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
            const TetMeshComponentDesc& desc);

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
            const TetMeshComponentDesc& desc);

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

        void RebuildTetFemCache();
        void EnsureFemSparsePattern(int worldNodeCount);
        void SyncCurrentPositionsFromWorld(const World& world);
        void SyncWorldTetsFromLocalTets(const World& world);

        int GetWorldNodeIndex(int localNodeIndex) const override;
        void SetLocalCurrentPosition(int localNodeIndex, const Vec3& position) override;

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
