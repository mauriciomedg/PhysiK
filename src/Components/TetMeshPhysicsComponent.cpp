#include "PhysiK/Components/TetMeshPhysicsComponent.h"

#include <algorithm>
#include <cmath>
#include <cstddef>

#include "PhysiK/Core/Solvers/SolverData.h"
#include "PhysiK/Core/World/World.h"

namespace PhysiK
{
    namespace
    {
        int FindLocalNodeIndex(
            const std::vector<int>& localToGlobalNodeIndex,
            int worldNodeIndex)
        {
            for (int localIndex = 0;
                 localIndex < static_cast<int>(localToGlobalNodeIndex.size());
                 ++localIndex)
            {
                if (localToGlobalNodeIndex[static_cast<std::size_t>(localIndex)] ==
                    worldNodeIndex)
                {
                    return localIndex;
                }
            }

            return -1;
        }

        Tet MakeTet(int node0, int node1, int node2, int node3)
        {
            Tet tet;
            tet.node0 = node0;
            tet.node1 = node1;
            tet.node2 = node2;
            tet.node3 = node3;
            return tet;
        }

        void ApplyMaterialToTet(Tet& tet, const Material& material)
        {
            tet.youngModulus = material.youngModulus;
            tet.poissonRatio = material.poissonRatio;
            tet.damping = material.damping;
        }

        int MapLocalToGlobal(
            const std::vector<int>& localToGlobalNodeIndex,
            int localNodeIndex)
        {
            if (localNodeIndex < 0 ||
                localNodeIndex >= static_cast<int>(localToGlobalNodeIndex.size()))
            {
                return -1;
            }

            return localToGlobalNodeIndex[static_cast<std::size_t>(localNodeIndex)];
        }

        std::vector<std::pair<int, int>> BuildGlobalSparsePattern(
            const std::vector<Tet>& localTets,
            const std::vector<int>& localToGlobalNodeIndex)
        {
            std::vector<std::pair<int, int>> globalCoordinates;
            const std::vector<std::pair<int, int>> localCoordinates =
                FEMModel::BuildSparsePatternFromTetConnectivity(localTets);
            globalCoordinates.reserve(localCoordinates.size());

            for (const std::pair<int, int>& coordinate : localCoordinates)
            {
                const int globalRow =
                    MapLocalToGlobal(localToGlobalNodeIndex, coordinate.first);
                const int globalColumn =
                    MapLocalToGlobal(localToGlobalNodeIndex, coordinate.second);
                if (globalRow < 0 || globalColumn < 0)
                {
                    continue;
                }

                globalCoordinates.push_back({globalRow, globalColumn});
            }

            return globalCoordinates;
        }

        void ApplyElementContribution(
            const TetElementContribution& contribution,
            SparseBlockMatrix& femSparseMatrix,
            SolverData& solverData,
            const std::vector<int>& localToGlobalNodeIndex)
        {
            // FEM data is local to the component; SolverData assembly must explicitly
            // map local node indices to global node indices.
            for (int node = 0; node < 4; ++node)
            {
                const int globalNodeIndex = MapLocalToGlobal(
                    localToGlobalNodeIndex,
                    contribution.localNodeIndices[node]);
                if (globalNodeIndex < 0)
                {
                    continue;
                }

                solverData.AddNodeForce(globalNodeIndex, contribution.forces[node]);
            }

            for (int rowNode = 0; rowNode < 4; ++rowNode)
            {
                const int globalRow = MapLocalToGlobal(
                    localToGlobalNodeIndex,
                    contribution.localNodeIndices[rowNode]);
                if (globalRow < 0)
                {
                    continue;
                }

                for (int columnNode = 0; columnNode < 4; ++columnNode)
                {
                    const int globalColumn = MapLocalToGlobal(
                        localToGlobalNodeIndex,
                        contribution.localNodeIndices[columnNode]);
                    if (globalColumn < 0)
                    {
                        continue;
                    }

                    femSparseMatrix.AddBlock(
                        globalRow,
                        globalColumn,
                        contribution.stiffness[rowNode][columnNode]);
                }
            }
        }

        void ApplyMassContribution(
            const TetMassContribution& contribution,
            const World& world,
            SolverData& solverData,
            const std::vector<int>& localToGlobalNodeIndex)
        {
            if (!std::isfinite(contribution.nodalMass) ||
                contribution.nodalMass <= 0.0f)
            {
                return;
            }

            for (int node = 0; node < 4; ++node)
            {
                const int globalNodeIndex = MapLocalToGlobal(
                    localToGlobalNodeIndex,
                    contribution.localNodeIndices[node]);
                if (globalNodeIndex < 0 ||
                    globalNodeIndex >= static_cast<int>(world.GetNodes().size()) ||
                    world.IsNodeFixed(globalNodeIndex))
                {
                    continue;
                }

                solverData.AddNodeMass(globalNodeIndex, contribution.nodalMass);
            }
        }
    }

    std::unique_ptr<TetMeshPhysicsComponent> TetMeshPhysicsComponent::CreateFromGlobalNodes(
        World& world,
        const int* globalNodeIndices,
        int nodeCount,
        const int* tetGlobalNodeIndices,
        int tetCount,
        const Material& material)
    {
        TetMeshPhysicsComponentDesc desc;
        desc.material = material;
        return CreateFromGlobalNodes(
            world,
            globalNodeIndices,
            nodeCount,
            tetGlobalNodeIndices,
            tetCount,
            desc);
    }

    std::unique_ptr<TetMeshPhysicsComponent> TetMeshPhysicsComponent::CreateFromGlobalNodes(
        World& world,
        const int* globalNodeIndices,
        int nodeCount,
        const int* tetGlobalNodeIndices,
        int tetCount,
        const TetMeshPhysicsComponentDesc& desc)
    {
        auto component = std::make_unique<TetMeshPhysicsComponent>();
        component->material = desc.material;
        component->selectedFemModel = desc.femModel;

        if (globalNodeIndices != nullptr && nodeCount > 0)
        {
            component->restNodePositions.reserve(static_cast<std::size_t>(nodeCount));
            component->nodePositions.reserve(static_cast<std::size_t>(nodeCount));
            component->localToGlobalNodeIndex.reserve(static_cast<std::size_t>(nodeCount));
            for (int i = 0; i < nodeCount; ++i)
            {
                const int worldNodeIndex = globalNodeIndices[i];
                if (worldNodeIndex < 0 ||
                    worldNodeIndex >= static_cast<int>(world.GetNodes().size()))
                {
                    continue;
                }

                component->localToGlobalNodeIndex.push_back(worldNodeIndex);
                component->restNodePositions.push_back(
                    world.GetNode(worldNodeIndex).restPosition);
                component->nodePositions.push_back(
                    world.GetNode(worldNodeIndex).position);
            }
        }

        if (tetGlobalNodeIndices != nullptr && tetCount > 0)
        {
            component->tets.reserve(static_cast<std::size_t>(tetCount));
            for (int i = 0; i < tetCount; ++i)
            {
                const int local0 = FindLocalNodeIndex(
                    component->localToGlobalNodeIndex,
                    tetGlobalNodeIndices[i * 4 + 0]);
                const int local1 = FindLocalNodeIndex(
                    component->localToGlobalNodeIndex,
                    tetGlobalNodeIndices[i * 4 + 1]);
                const int local2 = FindLocalNodeIndex(
                    component->localToGlobalNodeIndex,
                    tetGlobalNodeIndices[i * 4 + 2]);
                const int local3 = FindLocalNodeIndex(
                    component->localToGlobalNodeIndex,
                    tetGlobalNodeIndices[i * 4 + 3]);

                if (local0 < 0 || local1 < 0 || local2 < 0 || local3 < 0)
                {
                    continue;
                }

                component->tets.push_back(MakeTet(local0, local1, local2, local3));
            }
        }

        component->RebuildFemRestData();
        return component;
    }

    std::unique_ptr<TetMeshPhysicsComponent> TetMeshPhysicsComponent::CreateFromPositions(
        World& world,
        const Vec3* positions,
        const int* fixedNodeFlags,
        int nodeCount,
        const int* tetLocalNodeIndices,
        int tetCount,
        const Material& material)
    {
        TetMeshPhysicsComponentDesc desc;
        desc.material = material;
        return CreateFromPositions(
            world,
            positions,
            fixedNodeFlags,
            nodeCount,
            tetLocalNodeIndices,
            tetCount,
            desc);
    }

    std::unique_ptr<TetMeshPhysicsComponent> TetMeshPhysicsComponent::CreateFromPositions(
        World& world,
        const Vec3* positions,
        const int* fixedNodeFlags,
        int nodeCount,
        const int* tetLocalNodeIndices,
        int tetCount,
        const TetMeshPhysicsComponentDesc& desc)
    {
        auto component = std::make_unique<TetMeshPhysicsComponent>();
        component->material = desc.material;
        component->selectedFemModel = desc.femModel;

        if (positions != nullptr && nodeCount > 0)
        {
            component->localToGlobalNodeIndex.reserve(static_cast<std::size_t>(nodeCount));
            for (int i = 0; i < nodeCount; ++i)
            {
                const int nodeIndex = world.AddNode(positions[i]);
                if (fixedNodeFlags != nullptr && fixedNodeFlags[i] != 0)
                {
                    world.SetNodeFixed(nodeIndex, true);
                }
                component->localToGlobalNodeIndex.push_back(nodeIndex);
            }
        }

        component->SetGeometry(positions, nodeCount, tetLocalNodeIndices, tetCount);
        component->RebuildFemRestData();
        return component;
    }

    void TetMeshPhysicsComponent::SetMaterial(const Material& value)
    {
        material = value;
        for (Tet& tet : tets)
        {
            ApplyMaterialToTet(tet, material);
        }
        RebuildTetFemCache();
    }

    void TetMeshPhysicsComponent::RebuildFemRestData()
    {
        for (Tet& tet : tets)
        {
            ApplyMaterialToTet(tet, material);
            FEMModel::InitializeTetRestData(tet, restNodePositions);
        }

        RebuildTetFemCache();
        femSparsePatternDirty = true;
    }

    void TetMeshPhysicsComponent::RebuildTetFemCache()
    {
        tetFemCache.clear();
        tetFemCache.reserve(tets.size());
        for (const Tet& tet : tets)
        {
            tetFemCache.push_back(FEMModel::BuildTetFemCache(tet));
        }
    }

    void TetMeshPhysicsComponent::EnsureFemSparsePattern(int worldNodeCount)
    {
        if (!femSparsePatternDirty && femSparseMatrix.blockCount == worldNodeCount)
        {
            return;
        }

        femSparseMatrix.BuildPattern(
            worldNodeCount,
            BuildGlobalSparsePattern(tets, localToGlobalNodeIndex));
        femSparsePatternDirty = false;
    }

    void TetMeshPhysicsComponent::SyncCurrentPositionsFromWorld(const World& world)
    {
        nodePositions.resize(localToGlobalNodeIndex.size());
        nodeVelocities.resize(localToGlobalNodeIndex.size());
        for (int localIndex = 0;
             localIndex < static_cast<int>(localToGlobalNodeIndex.size());
             ++localIndex)
        {
            const int worldNodeIndex =
                localToGlobalNodeIndex[static_cast<std::size_t>(localIndex)];
            if (worldNodeIndex < 0 ||
                worldNodeIndex >= static_cast<int>(world.GetNodes().size()))
            {
                continue;
            }

            nodePositions[static_cast<std::size_t>(localIndex)] =
                world.GetNode(worldNodeIndex).position;
            nodeVelocities[static_cast<std::size_t>(localIndex)] =
                world.GetNode(worldNodeIndex).velocity;
        }
    }

    int TetMeshPhysicsComponent::GetGlobalNodeIndex(int localNodeIndex) const
    {
        if (localNodeIndex < 0 ||
            localNodeIndex >= static_cast<int>(localToGlobalNodeIndex.size()))
        {
            return -1;
        }

        return localToGlobalNodeIndex[static_cast<std::size_t>(localNodeIndex)];
    }

    void TetMeshPhysicsComponent::SetLocalCurrentPosition(
        int localNodeIndex,
        const Vec3& position)
    {
        TetMeshComponent::SetLocalCurrentPosition(localNodeIndex, position);
    }

    void TetMeshPhysicsComponent::SetTetActive(int tetIndex, bool active)
    {
        TetMeshComponent::SetTetActive(tetIndex, active);
    }

    void TetMeshPhysicsComponent::DeactivateTet(int tetIndex)
    {
        SetTetActive(tetIndex, false);
    }

    void TetMeshPhysicsComponent::UpdateSystem(
        World& world,
        SolverData& solverData,
        float dt)
    {
        (void)dt;
        SyncCurrentPositionsFromWorld(world);
        if (!physicsEnabled)
        {
            return;
        }

        std::vector<TetMassContribution> massContributions;
        FEMModel::ComputeLumpedMass(material, tets, massContributions);
        for (const TetMassContribution& contribution : massContributions)
        {
            ApplyMassContribution(
                contribution,
                world,
                solverData,
                localToGlobalNodeIndex);
        }

        EnsureFemSparsePattern(static_cast<int>(world.GetNodes().size()));
        if (tetFemCache.size() != tets.size())
        {
            RebuildTetFemCache();
        }

        std::vector<TetElementContribution> elementContributions;
        FEMModel::ComputeForces(
            GetFemModel(),
            tets,
            tetFemCache,
            nodePositions,
            nodeVelocities,
            elementContributions);

        femSparseMatrix.ClearValues();
        for (const TetElementContribution& contribution : elementContributions)
        {
            ApplyElementContribution(
                contribution,
                femSparseMatrix,
                solverData,
                localToGlobalNodeIndex);
        }

        for (int rowBlock = 0; rowBlock < femSparseMatrix.blockCount; ++rowBlock)
        {
            const int rowBegin = femSparseMatrix.rowStart[static_cast<std::size_t>(rowBlock)];
            const int rowEnd = femSparseMatrix.rowStart[static_cast<std::size_t>(rowBlock + 1)];
            for (int blockIndex = rowBegin; blockIndex < rowEnd; ++blockIndex)
            {
                solverData.AddStiffnessBlock(
                    rowBlock,
                    femSparseMatrix.colIndex[static_cast<std::size_t>(blockIndex)],
                    femSparseMatrix.values[static_cast<std::size_t>(blockIndex)]);
            }
        }
    }

    void TetMeshPhysicsComponent::PostUpdate(World& world, float dt)
    {
        if (dt > 0.0f)
        {
            SyncCurrentPositionsFromWorld(world);
        }
        TetMeshComponent::PostUpdate(world, dt);
    }
}
