#include "PhysiK/Components/TetMeshPhysicsComponent.h"

#include <algorithm>
#include <cstddef>

#include "PhysiK/Core/Solvers/SolverData.h"
#include "PhysiK/Core/World/World.h"

namespace PhysiK
{
    namespace
    {
        int FindLocalNodeIndex(const std::vector<int>& globalNodeIndices, int worldNodeIndex)
        {
            for (int localIndex = 0;
                 localIndex < static_cast<int>(globalNodeIndices.size());
                 ++localIndex)
            {
                if (globalNodeIndices[static_cast<std::size_t>(localIndex)] == worldNodeIndex)
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

        int MapLocalToWorld(
            const std::vector<int>& globalNodeIndices,
            int localNodeIndex)
        {
            if (localNodeIndex < 0)
            {
                return -1;
            }

            if (globalNodeIndices.empty())
            {
                return localNodeIndex;
            }

            if (localNodeIndex >= static_cast<int>(globalNodeIndices.size()))
            {
                return -1;
            }

            return globalNodeIndices[static_cast<std::size_t>(localNodeIndex)];
        }

        Tet BuildMappedTet(
            const Tet& localTet,
            const std::vector<int>& globalNodeIndices)
        {
            Tet worldTet = MakeTet(
                MapLocalToWorld(globalNodeIndices, localTet.node0),
                MapLocalToWorld(globalNodeIndices, localTet.node1),
                MapLocalToWorld(globalNodeIndices, localTet.node2),
                MapLocalToWorld(globalNodeIndices, localTet.node3));
            worldTet.active = localTet.active;
            return worldTet;
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
            component->globalNodeIndices.reserve(static_cast<std::size_t>(nodeCount));
            for (int i = 0; i < nodeCount; ++i)
            {
                const int worldNodeIndex = globalNodeIndices[i];
                if (worldNodeIndex < 0 ||
                    worldNodeIndex >= static_cast<int>(world.GetNodes().size()))
                {
                    continue;
                }

                component->globalNodeIndices.push_back(worldNodeIndex);
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
                    component->globalNodeIndices,
                    tetGlobalNodeIndices[i * 4 + 0]);
                const int local1 = FindLocalNodeIndex(
                    component->globalNodeIndices,
                    tetGlobalNodeIndices[i * 4 + 1]);
                const int local2 = FindLocalNodeIndex(
                    component->globalNodeIndices,
                    tetGlobalNodeIndices[i * 4 + 2]);
                const int local3 = FindLocalNodeIndex(
                    component->globalNodeIndices,
                    tetGlobalNodeIndices[i * 4 + 3]);

                if (local0 < 0 || local1 < 0 || local2 < 0 || local3 < 0)
                {
                    continue;
                }

                component->tets.push_back(MakeTet(local0, local1, local2, local3));
            }
        }

        component->RebuildWorldTets(world);
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
            component->globalNodeIndices.reserve(static_cast<std::size_t>(nodeCount));
            for (int i = 0; i < nodeCount; ++i)
            {
                const int nodeIndex = world.AddNode(positions[i]);
                if (fixedNodeFlags != nullptr && fixedNodeFlags[i] != 0)
                {
                    world.SetNodeFixed(nodeIndex, true);
                }
                component->globalNodeIndices.push_back(nodeIndex);
            }
        }

        component->SetGeometry(positions, nodeCount, tetLocalNodeIndices, tetCount);
        component->RebuildWorldTets(world);
        return component;
    }

    void TetMeshPhysicsComponent::SetMaterial(const Material& value)
    {
        material = value;
        for (Tet& tet : globalTets)
        {
            ApplyMaterialToTet(tet, material);
        }
        RebuildTetFemCache();
    }

    void TetMeshPhysicsComponent::RebuildWorldTets(const World& world)
    {
        globalTets.clear();
        globalTets.reserve(tets.size());
        for (const Tet& localTet : tets)
        {
            Tet mappedTet = BuildMappedTet(localTet, globalNodeIndices);
            ApplyMaterialToTet(mappedTet, material);
            FEMModel::InitializeTetRestData(mappedTet, world.GetNodes());
            globalTets.push_back(mappedTet);
        }

        RebuildTetFemCache();
        femSparsePatternDirty = true;
    }

    void TetMeshPhysicsComponent::RebuildTetFemCache()
    {
        tetFemCache.clear();
        tetFemCache.reserve(globalTets.size());
        for (const Tet& tet : globalTets)
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
            FEMModel::BuildSparsePatternFromTetConnectivity(globalTets));
        femSparsePatternDirty = false;
    }

    void TetMeshPhysicsComponent::SyncWorldTetActiveStates()
    {
        const std::size_t count = std::min(tets.size(), globalTets.size());
        for (std::size_t tetIndex = 0; tetIndex < count; ++tetIndex)
        {
            globalTets[tetIndex].active = tets[tetIndex].active;
        }
    }

    void TetMeshPhysicsComponent::SyncCurrentPositionsFromWorld(const World& world)
    {
        nodePositions.resize(globalNodeIndices.size());
        for (int localIndex = 0;
             localIndex < static_cast<int>(globalNodeIndices.size());
             ++localIndex)
        {
            const int worldNodeIndex = globalNodeIndices[static_cast<std::size_t>(localIndex)];
            if (worldNodeIndex < 0 ||
                worldNodeIndex >= static_cast<int>(world.GetNodes().size()))
            {
                continue;
            }

            nodePositions[static_cast<std::size_t>(localIndex)] =
                world.GetNode(worldNodeIndex).position;
        }
    }

    int TetMeshPhysicsComponent::GetGlobalNodeIndex(int localNodeIndex) const
    {
        if (localNodeIndex < 0 ||
            localNodeIndex >= static_cast<int>(globalNodeIndices.size()))
        {
            return -1;
        }

        return globalNodeIndices[static_cast<std::size_t>(localNodeIndex)];
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
        if (tetIndex < 0 || tetIndex >= static_cast<int>(globalTets.size()))
        {
            return;
        }

        globalTets[static_cast<std::size_t>(tetIndex)].active = active;
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

        SyncWorldTetActiveStates();
        FEMModel::AssembleLumpedMass(material, globalTets, world, solverData);
        EnsureFemSparsePattern(static_cast<int>(world.GetNodes().size()));
        if (tetFemCache.size() != globalTets.size())
        {
            RebuildTetFemCache();
        }

        SolverData femSolverData;
        FEMModel::AccumulateForces(
            GetFemModel(),
            globalTets,
            tetFemCache,
            world.GetNodes(),
            femSolverData);

        for (const SolverData::NodeForce& force : femSolverData.GetNodeForces())
        {
            solverData.AddNodeForce(force.node, force.force);
        }

        femSparseMatrix.ClearValues();
        for (const SolverData::StiffnessBlock& block : femSolverData.GetStiffnessBlocks())
        {
            femSparseMatrix.AddBlock(block.nodeA, block.nodeB, block.block);
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
        SyncCurrentPositionsFromWorld(world);
        TetMeshComponent::PostUpdate(world, dt);
    }
}
