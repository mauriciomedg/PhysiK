#include "PhysiK/Components/TetMeshComponent.h"

#include <algorithm>
#include <cmath>
#include <optional>

#include "PhysiK/Core/Solvers/SolverData.h"
#include "PhysiK/Core/World/World.h"
#if defined(PHYSIK_ENABLE_PERF_LOGGING)
#include "PhysiK/Core/Performance/PerformanceLogger.h"
#endif

namespace PhysiK
{
    namespace
    {
        void AppendTetFromGlobalNodes(
            TetMeshComponent& component,
            World& world,
            int node0,
            int node1,
            int node2,
            int node3)
        {
            Tet tet;
            tet.node0 = node0;
            tet.node1 = node1;
            tet.node2 = node2;
            tet.node3 = node3;
            tet.youngModulus = component.material.youngModulus;
            tet.poissonRatio = component.material.poissonRatio;
            tet.damping = component.material.damping;
            FEMModel::InitializeTetRestData(tet, world.GetNodes());
            component.tets.push_back(tet);
            component.tetFemCache.push_back(FEMModel::BuildTetFemCache(component.tets.back()));
            component.MarkFemSparsePatternDirty();
        }

        void AddLumpedMassToSolverData(
            const World& world,
            SolverData& solverData,
            int nodeIndex,
            float mass)
        {
            if (nodeIndex < 0 || !std::isfinite(mass))
            {
                return;
            }

            if (world.IsNodeFixed(nodeIndex))
            {
                return;
            }

            solverData.AddNodeMass(nodeIndex, std::max(0.0f, mass));
        }

        void AssembleLumpedMass(
            const TetMeshComponent& component,
            const World& world,
            SolverData& solverData)
        {
            const float density = std::max(0.0f, component.material.density);
            if (!std::isfinite(density))
            {
                return;
            }

            for (const Tet& tet : component.tets)
            {
                if (!tet.active)
                {
                    continue;
                }

                if (!std::isfinite(tet.restVolume) || tet.restVolume <= 0.0f)
                {
                    continue;
                }

                const float nodalMass = density * tet.restVolume * 0.25f;
                AddLumpedMassToSolverData(world, solverData, tet.node0, nodalMass);
                AddLumpedMassToSolverData(world, solverData, tet.node1, nodalMass);
                AddLumpedMassToSolverData(world, solverData, tet.node2, nodalMass);
                AddLumpedMassToSolverData(world, solverData, tet.node3, nodalMass);
            }
        }

        std::vector<std::pair<int, int>> BuildSparsePatternFromTetConnectivity(
            const std::vector<Tet>& tets)
        {
            std::vector<std::pair<int, int>> blockCoordinates;
            blockCoordinates.reserve(tets.size() * 16u);

            for (const Tet& tet : tets)
            {
                const int nodes[4] = {tet.node0, tet.node1, tet.node2, tet.node3};
                for (int row = 0; row < 4; ++row)
                {
                    for (int column = 0; column < 4; ++column)
                    {
                        blockCoordinates.push_back({nodes[row], nodes[column]});
                    }
                }
            }

            return blockCoordinates;
        }
    }

    std::unique_ptr<TetMeshComponent> TetMeshComponent::CreateFromGlobalNodes(
        World& world,
        const int* globalNodeIndices,
        int nodeCount,
        const int* tetGlobalNodeIndices,
        int tetCount,
        const Material& material)
    {
        TetMeshComponentDesc desc;
        desc.material = material;
        return CreateFromGlobalNodes(
            world,
            globalNodeIndices,
            nodeCount,
            tetGlobalNodeIndices,
            tetCount,
            desc);
    }

    std::unique_ptr<TetMeshComponent> TetMeshComponent::CreateFromGlobalNodes(
        World& world,
        const int* globalNodeIndices,
        int nodeCount,
        const int* tetGlobalNodeIndices,
        int tetCount,
        const TetMeshComponentDesc& desc)
    {
        auto component = std::make_unique<TetMeshComponent>();
        component->material = desc.material;
        component->selectedFemModel = desc.femModel;

        if (globalNodeIndices != nullptr && nodeCount > 0)
        {
            component->nodeIndices.assign(globalNodeIndices, globalNodeIndices + nodeCount);
        }

        if (tetGlobalNodeIndices != nullptr && tetCount > 0)
        {
            component->tets.reserve(static_cast<std::size_t>(tetCount));
            component->tetFemCache.reserve(static_cast<std::size_t>(tetCount));
            for (int i = 0; i < tetCount; ++i)
            {
                AppendTetFromGlobalNodes(
                    *component,
                    world,
                    tetGlobalNodeIndices[i * 4 + 0],
                    tetGlobalNodeIndices[i * 4 + 1],
                    tetGlobalNodeIndices[i * 4 + 2],
                    tetGlobalNodeIndices[i * 4 + 3]);
            }
        }

        return component;
    }

    std::unique_ptr<TetMeshComponent> TetMeshComponent::CreateFromPositions(
        World& world,
        const Vec3* positions,
        const int* fixedNodeFlags,
        int nodeCount,
        const int* tetLocalNodeIndices,
        int tetCount,
        const Material& material)
    {
        TetMeshComponentDesc desc;
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

    std::unique_ptr<TetMeshComponent> TetMeshComponent::CreateFromPositions(
        World& world,
        const Vec3* positions,
        const int* fixedNodeFlags,
        int nodeCount,
        const int* tetLocalNodeIndices,
        int tetCount,
        const TetMeshComponentDesc& desc)
    {
        auto component = std::make_unique<TetMeshComponent>();
        component->material = desc.material;
        component->selectedFemModel = desc.femModel;

        if (positions != nullptr && nodeCount > 0)
        {
            component->nodeIndices.reserve(static_cast<std::size_t>(nodeCount));
            for (int i = 0; i < nodeCount; ++i)
            {
                const int nodeIndex = world.AddNode(positions[i]);
                if (fixedNodeFlags != nullptr && fixedNodeFlags[i] != 0)
                {
                    world.SetNodeFixed(nodeIndex, true);
                }
                component->nodeIndices.push_back(nodeIndex);
            }
        }

        if (tetLocalNodeIndices != nullptr && tetCount > 0)
        {
            component->tets.reserve(static_cast<std::size_t>(tetCount));
            component->tetFemCache.reserve(static_cast<std::size_t>(tetCount));
            for (int i = 0; i < tetCount; ++i)
            {
                const int local0 = tetLocalNodeIndices[i * 4 + 0];
                const int local1 = tetLocalNodeIndices[i * 4 + 1];
                const int local2 = tetLocalNodeIndices[i * 4 + 2];
                const int local3 = tetLocalNodeIndices[i * 4 + 3];
                const int nodeCountInComponent = static_cast<int>(component->nodeIndices.size());

                if (local0 < 0 || local0 >= nodeCountInComponent ||
                    local1 < 0 || local1 >= nodeCountInComponent ||
                    local2 < 0 || local2 >= nodeCountInComponent ||
                    local3 < 0 || local3 >= nodeCountInComponent)
                {
                    continue;
                }

                AppendTetFromGlobalNodes(
                    *component,
                    world,
                    component->nodeIndices[static_cast<std::size_t>(local0)],
                    component->nodeIndices[static_cast<std::size_t>(local1)],
                    component->nodeIndices[static_cast<std::size_t>(local2)],
                    component->nodeIndices[static_cast<std::size_t>(local3)]);
            }
        }

        return component;
    }

    void TetMeshComponent::SetMaterial(const Material& value)
    {
        material = value;
        for (Tet& tet : tets)
        {
            tet.youngModulus = material.youngModulus;
            tet.poissonRatio = material.poissonRatio;
            tet.damping = material.damping;
        }
        RebuildTetFemCache();
    }

    void TetMeshComponent::RebuildTetFemCache()
    {
        tetFemCache.clear();
        tetFemCache.reserve(tets.size());
        for (const Tet& tet : tets)
        {
            tetFemCache.push_back(FEMModel::BuildTetFemCache(tet));
        }
    }

    void TetMeshComponent::EnsureFemSparsePattern(int worldNodeCount)
    {
        if (!femSparsePatternDirty && femSparseMatrix.blockCount == worldNodeCount)
        {
            return;
        }

        femSparseMatrix.BuildPattern(
            worldNodeCount,
            BuildSparsePatternFromTetConnectivity(tets));
        femSparsePatternDirty = false;
    }

    bool TetMeshComponent::IsTetActive(int tetIndex) const
    {
        if (tetIndex < 0 || tetIndex >= static_cast<int>(tets.size()))
        {
            return false;
        }

        return tets[static_cast<std::size_t>(tetIndex)].active;
    }

    void TetMeshComponent::SetTetActive(int tetIndex, bool active)
    {
        if (tetIndex < 0 || tetIndex >= static_cast<int>(tets.size()))
        {
            return;
        }

        tets[static_cast<std::size_t>(tetIndex)].active = active;
    }

    void TetMeshComponent::DeactivateTet(int tetIndex)
    {
        SetTetActive(tetIndex, false);
        topologyDirty = true;
    }

    int TetMeshComponent::GetActiveTetCount() const
    {
        int activeCount = 0;
        for (const Tet& tet : tets)
        {
            if (tet.active)
            {
                ++activeCount;
            }
        }

        return activeCount;
    }

    void TetMeshComponent::UpdateSystem(
        World& world,
        SolverData& solverData,
        float dt)
    {
        AssembleLumpedMass(*this, world, solverData);
        EnsureFemSparsePattern(static_cast<int>(world.GetNodes().size()));
        if (tetFemCache.size() != tets.size())
        {
            RebuildTetFemCache();
        }

        SolverData femSolverData;
        femModel.UpdateSystem(world, *this, femSolverData, dt);

        for (const SolverData::NodeForce& force : femSolverData.GetNodeForces())
        {
            solverData.AddNodeForce(force.node, force.force);
        }

        femSparseMatrix.ClearValues();
#if defined(PHYSIK_ENABLE_PERF_LOGGING)
        PerformanceLogRecord* performanceRecord = FEMModel::GetPerformanceLogRecord();
        const bool logPerformance = performanceRecord != nullptr;
        std::optional<PerformanceTimer> matrixAddBlockTimer;
        if (logPerformance)
        {
            matrixAddBlockTimer.emplace();
        }
#endif
        for (const SolverData::StiffnessBlock& block : femSolverData.GetStiffnessBlocks())
        {
            femSparseMatrix.AddBlock(block.nodeA, block.nodeB, block.block);
        }
#if defined(PHYSIK_ENABLE_PERF_LOGGING)
        if (logPerformance)
        {
            const double matrixAddBlockMilliseconds =
                matrixAddBlockTimer->ElapsedMilliseconds();
            performanceRecord->assembleMatrixAddBlockMs += matrixAddBlockMilliseconds;
            performanceRecord->tetMatrixWriteMs += matrixAddBlockMilliseconds;
        }
#endif

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

    void TetMeshComponent::Execute(World& world)
    {
        if (!topologyDirty)
            return;

        PhysicsEvent event;
        event.type = PhysicsEventType::TetMeshTopologyChanged;
        event.world = &world;
        event.sender = this;
        //event.objectIndex = lastChangedTetIndex; // optional

        world.EmitEvent(event);

        topologyDirty = false;
    }
}
