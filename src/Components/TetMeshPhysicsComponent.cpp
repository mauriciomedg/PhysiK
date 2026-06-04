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
                *this);
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
                *this);
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
