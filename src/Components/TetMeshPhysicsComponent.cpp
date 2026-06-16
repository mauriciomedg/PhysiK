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
        void ApplyMaterialToTet(Tet& tet, const Material& material)
        {
            tet.youngModulus = material.youngModulus;
            tet.poissonRatio = material.poissonRatio;
            tet.damping = material.damping;
        }

        void AddGeneratedNodesToWorld(
            World& world,
            TetMeshPhysicsComponent& component)
        {
            component.globalNodeBeginIndex =
                static_cast<int>(world.GetNodes().size());
            component.globalNodeCount =
                static_cast<int>(component.restNodePositions.size());

            for (const Vec3& position : component.restNodePositions)
            {
                world.AddNode(position);
            }
        }

        std::vector<std::pair<int, int>> BuildGlobalSparsePattern(
            const std::vector<Tet>& localTets,
            const TetMeshPhysicsComponent& component)
        {
            std::vector<std::pair<int, int>> globalCoordinates;
            const std::vector<std::pair<int, int>> localCoordinates =
                FEMModel::BuildSparsePatternFromTetConnectivity(localTets);
            globalCoordinates.reserve(localCoordinates.size());

            for (const std::pair<int, int>& coordinate : localCoordinates)
            {
                const int globalRow =
                    component.GetGlobalNodeIndex(coordinate.first);
                const int globalColumn =
                    component.GetGlobalNodeIndex(coordinate.second);
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
            const TetMeshPhysicsComponent& component)
        {
            // FEM data is local to the component; SolverData assembly must explicitly
            // map local node indices to global node indices.
            for (int node = 0; node < 4; ++node)
            {
                const int globalNodeIndex = component.GetGlobalNodeIndex(
                    contribution.localNodeIndices[node]);
                if (globalNodeIndex < 0)
                {
                    continue;
                }

                solverData.AddNodeForce(globalNodeIndex, contribution.forces[node]);
            }

            for (int rowNode = 0; rowNode < 4; ++rowNode)
            {
                const int globalRow = component.GetGlobalNodeIndex(
                    contribution.localNodeIndices[rowNode]);
                if (globalRow < 0)
                {
                    continue;
                }

                for (int columnNode = 0; columnNode < 4; ++columnNode)
                {
                    const int globalColumn = component.GetGlobalNodeIndex(
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

        void ApplyElementContribution(
            const TetElementContribution& contribution,
            SparseBlockMatrix& femSparseMatrix,
            SolverData& solverData,
            const std::vector<int>& cachedLocalToGlobalNodeIndices)
        {
            for (int node = 0; node < 4; ++node)
            {
                const int localNodeIndex = contribution.localNodeIndices[node];
                if (localNodeIndex < 0 ||
                    localNodeIndex >= static_cast<int>(cachedLocalToGlobalNodeIndices.size()))
                {
                    continue;
                }

                const int globalNodeIndex =
                    cachedLocalToGlobalNodeIndices[static_cast<std::size_t>(localNodeIndex)];
                if (globalNodeIndex < 0)
                {
                    continue;
                }

                solverData.AddNodeForce(globalNodeIndex, contribution.forces[node]);
            }

            for (int rowNode = 0; rowNode < 4; ++rowNode)
            {
                const int localRow = contribution.localNodeIndices[rowNode];
                if (localRow < 0 ||
                    localRow >= static_cast<int>(cachedLocalToGlobalNodeIndices.size()))
                {
                    continue;
                }

                const int globalRow =
                    cachedLocalToGlobalNodeIndices[static_cast<std::size_t>(localRow)];
                if (globalRow < 0)
                {
                    continue;
                }

                for (int columnNode = 0; columnNode < 4; ++columnNode)
                {
                    const int localColumn = contribution.localNodeIndices[columnNode];
                    if (localColumn < 0 ||
                        localColumn >= static_cast<int>(cachedLocalToGlobalNodeIndices.size()))
                    {
                        continue;
                    }

                    const int globalColumn =
                        cachedLocalToGlobalNodeIndices[static_cast<std::size_t>(localColumn)];
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
            const TetMeshPhysicsComponent& component)
        {
            if (!std::isfinite(contribution.nodalMass) ||
                contribution.nodalMass <= 0.0f)
            {
                return;
            }

            for (int node = 0; node < 4; ++node)
            {
                const int globalNodeIndex = component.GetGlobalNodeIndex(
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

    ComponentExecutionPriority
    TetMeshPhysicsComponent::GetExecutionPriority() const
    {
        return ComponentExecutionPriority::
            TetMeshPhysicsComponent;
    }

    std::unique_ptr<TetMeshPhysicsComponent>
    TetMeshPhysicsComponent::CreateFromGeneratedTetMesh(
        World& world,
        const GeneratedTetMesh& generatedMesh,
        const Material& material)
    {
        TetMeshPhysicsComponentDesc desc;
        desc.material = material;
        return CreateFromGeneratedTetMesh(world, generatedMesh, desc);
    }

    std::unique_ptr<TetMeshPhysicsComponent>
    TetMeshPhysicsComponent::CreateFromGeneratedTetMesh(
        World& world,
        const GeneratedTetMesh& generatedMesh,
        const TetMeshPhysicsComponentDesc& desc)
    {
        auto component = std::make_unique<TetMeshPhysicsComponent>();
        component->material = desc.material;
        component->selectedFemModel = desc.femModel;

        component->SetGeometry(generatedMesh);
        AddGeneratedNodesToWorld(world, *component);
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
        MarkFEMCacheDirty();
    }

    void TetMeshPhysicsComponent::MarkFEMCacheDirty()
    {
        femCacheDirty = true;
    }

    void TetMeshPhysicsComponent::RebuildFEMCacheIfNeeded(const World& world)
    {
        if (!femCacheDirty)
        {
#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
            ++femCacheReuseCount;
#endif
            return;
        }

        RebuildFEMCache(world);
    }

    void TetMeshPhysicsComponent::RebuildFEMCache(const World& world)
    {
        (void)world;

        cachedTetEntries.clear();
        cachedActiveTets.clear();
        cachedActiveTetFemCache.clear();
        cachedLocalToGlobalNodeIndices.assign(
            static_cast<std::size_t>(std::max(0, globalNodeCount)),
            -1);

        for (int localNodeIndex = 0; localNodeIndex < globalNodeCount; ++localNodeIndex)
        {
            cachedLocalToGlobalNodeIndices[static_cast<std::size_t>(localNodeIndex)] =
                GetGlobalNodeIndex(localNodeIndex);
        }

        const float density = std::max(0.0f, material.density);
        if (!std::isfinite(density))
        {
            femCacheDirty = false;
#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
            ++femCacheRebuildCount;
#endif
            return;
        }

        cachedTetEntries.reserve(tets.size());
        cachedActiveTets.reserve(tets.size());
        cachedActiveTetFemCache.reserve(tets.size());

        for (int tetIndex = 0; tetIndex < static_cast<int>(tets.size()); ++tetIndex)
        {
            const Tet& tet = tets[static_cast<std::size_t>(tetIndex)];
            if (!tet.active)
            {
                continue;
            }

            const int localNodeIndices[4] = {tet.node0, tet.node1, tet.node2, tet.node3};
            CachedTetEntry entry;
            entry.tetIndex = tetIndex;
            bool hasValidGlobalNodes = true;
            for (int node = 0; node < 4; ++node)
            {
                entry.localNodeIndices[node] = localNodeIndices[node];
                entry.globalNodeIndices[node] = GetCachedGlobalNodeIndex(localNodeIndices[node]);
                hasValidGlobalNodes = hasValidGlobalNodes && entry.globalNodeIndices[node] >= 0;
            }

            for (int row = 0; row < 4; ++row)
            {
                for (int column = 0; column < 4; ++column)
                {
                    entry.stiffnessNodePairs[row * 4 + column] = {
                        entry.globalNodeIndices[row],
                        entry.globalNodeIndices[column]};
                }
            }

            if (!hasValidGlobalNodes ||
                !std::isfinite(tet.restVolume) ||
                tet.restVolume <= 0.0f)
            {
                continue;
            }

            entry.nodalMass = density * tet.restVolume * 0.25f;
            entry.femCache = FEMModel::BuildTetFemCache(tet);

            cachedTetEntries.push_back(entry);
            cachedActiveTets.push_back(tet);
            cachedActiveTetFemCache.push_back(entry.femCache);
        }

        femCacheDirty = false;
#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
        ++femCacheRebuildCount;
#endif
    }

    void TetMeshPhysicsComponent::RebuildFemRestData()
    {
        for (Tet& tet : tets)
        {
            ApplyMaterialToTet(tet, material);
            FEMModel::InitializeTetRestData(tet, restNodePositions);
        }

        RebuildTetFemCache();
        MarkFEMCacheDirty();
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
        MarkFEMCacheDirty();
    }

    void TetMeshPhysicsComponent::EnsureFemSparsePattern(int worldNodeCount)
    {
        if (!femSparsePatternDirty && femSparseMatrix.blockCount == worldNodeCount)
        {
            return;
        }

        femSparseMatrix.BuildPattern(
            worldNodeCount,
            BuildGlobalSparsePattern(tets, *this));
        femSparsePatternDirty = false;
    }

    void TetMeshPhysicsComponent::SyncCurrentPositionsFromWorld(const World& world)
    {
        nodePositions.resize(static_cast<std::size_t>(globalNodeCount));
        nodeVelocities.resize(static_cast<std::size_t>(globalNodeCount));
        for (int localIndex = 0; localIndex < globalNodeCount; ++localIndex)
        {
            const int worldNodeIndex = GetGlobalNodeIndex(localIndex);
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

    int TetMeshPhysicsComponent::GetGlobalNodeBeginIndex() const
    {
        return globalNodeBeginIndex;
    }

    int TetMeshPhysicsComponent::GetGlobalNodeCount() const
    {
        return globalNodeCount;
    }

    int TetMeshPhysicsComponent::GetGlobalNodeIndex(int localNodeIndex) const
    {
        if (globalNodeBeginIndex < 0 ||
            localNodeIndex < 0 ||
            localNodeIndex >= globalNodeCount)
        {
            return -1;
        }

        return globalNodeBeginIndex + localNodeIndex;
    }

    void TetMeshPhysicsComponent::SetLocalCurrentPosition(
        int localNodeIndex,
        const Vec3& position)
    {
        TetMeshComponent::SetLocalCurrentPosition(localNodeIndex, position);
    }

    bool TetMeshPhysicsComponent::SetTetActive(int tetIndex, bool active)
    {
        const bool changed = TetMeshComponent::SetTetActive(tetIndex, active);
        if (changed)
        {
            MarkFEMCacheDirty();
        }
        return changed;
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

        RebuildFEMCacheIfNeeded(world);

        for (const CachedTetEntry& entry : cachedTetEntries)
        {
            if (!std::isfinite(entry.nodalMass) || entry.nodalMass <= 0.0f)
            {
                continue;
            }

            for (int node = 0; node < 4; ++node)
            {
                const int globalNodeIndex = entry.globalNodeIndices[node];
                if (globalNodeIndex < 0 ||
                    globalNodeIndex >= static_cast<int>(world.GetNodes().size()) ||
                    world.IsNodeFixed(globalNodeIndex))
                {
                    continue;
                }

                solverData.AddNodeMass(globalNodeIndex, entry.nodalMass);
            }
        }

        EnsureFemSparsePattern(static_cast<int>(world.GetNodes().size()));

        std::vector<TetElementContribution> elementContributions;
        FEMModel::ComputeForces(
            GetFemModel(),
            cachedActiveTets,
            cachedActiveTetFemCache,
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
                cachedLocalToGlobalNodeIndices);
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

    int TetMeshPhysicsComponent::GetCachedGlobalNodeIndex(int localNodeIndex) const
    {
        if (localNodeIndex < 0 ||
            localNodeIndex >= static_cast<int>(cachedLocalToGlobalNodeIndices.size()))
        {
            return -1;
        }

        return cachedLocalToGlobalNodeIndices[static_cast<std::size_t>(localNodeIndex)];
    }
}
