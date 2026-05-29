#include "PhysiK/Components/TetMeshPhysicsComponent.h"

#include <algorithm>
#include <cmath>
#include <cstddef>

#include "PhysiK/Components/TetMeshPreprocessor.h"
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

        std::vector<int> rawLocalToGlobalNodeIndex;
        std::vector<Vec3> rawPositions;
        if (globalNodeIndices != nullptr && nodeCount > 0)
        {
            rawPositions.reserve(static_cast<std::size_t>(nodeCount));
            rawLocalToGlobalNodeIndex.reserve(static_cast<std::size_t>(nodeCount));
            for (int i = 0; i < nodeCount; ++i)
            {
                const int worldNodeIndex = globalNodeIndices[i];
                if (worldNodeIndex < 0 ||
                    worldNodeIndex >= static_cast<int>(world.GetNodes().size()))
                {
                    continue;
                }

                rawLocalToGlobalNodeIndex.push_back(worldNodeIndex);
                rawPositions.push_back(world.GetNode(worldNodeIndex).restPosition);
            }
        }

        std::vector<int> rawTetLocalNodeIndices;
        if (tetGlobalNodeIndices != nullptr && tetCount > 0)
        {
            rawTetLocalNodeIndices.reserve(static_cast<std::size_t>(tetCount) * 4u);
            for (int i = 0; i < tetCount; ++i)
            {
                const int local0 = FindLocalNodeIndex(
                    rawLocalToGlobalNodeIndex,
                    tetGlobalNodeIndices[i * 4 + 0]);
                const int local1 = FindLocalNodeIndex(
                    rawLocalToGlobalNodeIndex,
                    tetGlobalNodeIndices[i * 4 + 1]);
                const int local2 = FindLocalNodeIndex(
                    rawLocalToGlobalNodeIndex,
                    tetGlobalNodeIndices[i * 4 + 2]);
                const int local3 = FindLocalNodeIndex(
                    rawLocalToGlobalNodeIndex,
                    tetGlobalNodeIndices[i * 4 + 3]);

                if (local0 < 0 || local1 < 0 || local2 < 0 || local3 < 0)
                {
                    continue;
                }

                rawTetLocalNodeIndices.push_back(local0);
                rawTetLocalNodeIndices.push_back(local1);
                rawTetLocalNodeIndices.push_back(local2);
                rawTetLocalNodeIndices.push_back(local3);
            }
        }

        const TetMeshPreprocessResult preprocessed =
            PreprocessTetMesh(
                rawPositions.data(),
                static_cast<int>(rawPositions.size()),
                rawTetLocalNodeIndices.data(),
                static_cast<int>(rawTetLocalNodeIndices.size() / 4u));
        component->restNodePositions = preprocessed.positions;
        component->nodePositions = component->restNodePositions;
        component->localToGlobalNodeIndex.resize(preprocessed.positions.size(), -1);
        for (int newNode = 0;
             newNode < static_cast<int>(preprocessed.newNodeToFirstOldNode.size());
             ++newNode)
        {
            const int oldNode =
                preprocessed.newNodeToFirstOldNode[static_cast<std::size_t>(newNode)];
            if (oldNode >= 0 &&
                oldNode < static_cast<int>(rawLocalToGlobalNodeIndex.size()))
            {
                component->localToGlobalNodeIndex[static_cast<std::size_t>(newNode)] =
                    rawLocalToGlobalNodeIndex[static_cast<std::size_t>(oldNode)];
            }
        }

        const int finalTetCount =
            static_cast<int>(preprocessed.tetLocalNodeIndices.size() / 4u);
        component->nodePositions.resize(component->localToGlobalNodeIndex.size());
        for (int localNode = 0;
             localNode < static_cast<int>(component->localToGlobalNodeIndex.size());
             ++localNode)
        {
            const int worldNode =
                component->localToGlobalNodeIndex[static_cast<std::size_t>(localNode)];
            if (worldNode >= 0 &&
                worldNode < static_cast<int>(world.GetNodes().size()))
            {
                component->nodePositions[static_cast<std::size_t>(localNode)] =
                    world.GetNode(worldNode).position;
            }
        }

        component->tets.reserve(static_cast<std::size_t>(finalTetCount));
        for (int tetIndex = 0; tetIndex < finalTetCount; ++tetIndex)
        {
            component->tets.push_back(MakeTet(
                preprocessed.tetLocalNodeIndices[tetIndex * 4 + 0],
                preprocessed.tetLocalNodeIndices[tetIndex * 4 + 1],
                preprocessed.tetLocalNodeIndices[tetIndex * 4 + 2],
                preprocessed.tetLocalNodeIndices[tetIndex * 4 + 3]));
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

        const TetMeshPreprocessResult preprocessed =
            PreprocessTetMesh(positions, nodeCount, tetLocalNodeIndices, tetCount);

        if (!preprocessed.positions.empty())
        {
            component->localToGlobalNodeIndex.reserve(preprocessed.positions.size());
            for (int newNode = 0;
                 newNode < static_cast<int>(preprocessed.positions.size());
                 ++newNode)
            {
                const int nodeIndex = world.AddNode(
                    preprocessed.positions[static_cast<std::size_t>(newNode)]);

                bool fixed = false;
                if (fixedNodeFlags != nullptr)
                {
                    for (int oldNode = 0; oldNode < nodeCount; ++oldNode)
                    {
                        if (preprocessed.oldNodeToNewNode[static_cast<std::size_t>(oldNode)] ==
                                newNode &&
                            fixedNodeFlags[oldNode] != 0)
                        {
                            fixed = true;
                            break;
                        }
                    }
                }

                if (fixed)
                {
                    world.SetNodeFixed(nodeIndex, true);
                }

                component->localToGlobalNodeIndex.push_back(nodeIndex);
            }
        }

        component->restNodePositions = preprocessed.positions;
        component->nodePositions = component->restNodePositions;
        const int finalTetCount =
            static_cast<int>(preprocessed.tetLocalNodeIndices.size() / 4u);
        component->tets.reserve(static_cast<std::size_t>(finalTetCount));
        for (int tetIndex = 0; tetIndex < finalTetCount; ++tetIndex)
        {
            component->tets.push_back(MakeTet(
                preprocessed.tetLocalNodeIndices[tetIndex * 4 + 0],
                preprocessed.tetLocalNodeIndices[tetIndex * 4 + 1],
                preprocessed.tetLocalNodeIndices[tetIndex * 4 + 2],
                preprocessed.tetLocalNodeIndices[tetIndex * 4 + 3]));
        }
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

    bool TetMeshPhysicsComponent::SetTetActive(int tetIndex, bool active)
    {
        return TetMeshComponent::SetTetActive(tetIndex, active);
    }

    bool TetMeshPhysicsComponent::DeactivateTet(int tetIndex)
    {
        return SetTetActive(tetIndex, false);
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
