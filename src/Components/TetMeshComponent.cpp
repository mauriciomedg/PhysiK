#include "PhysiK/Components/TetMeshComponent.h"

#include <algorithm>
#include <cstddef>

#include "PhysiK/Core/Solvers/SolverData.h"
#include "PhysiK/Core/World/World.h"

namespace PhysiK
{
    namespace
    {
        const Vec3 ZeroVec3;

        int FindLocalNodeIndex(const std::vector<int>& worldNodeIndices, int worldNodeIndex)
        {
            for (int localIndex = 0;
                 localIndex < static_cast<int>(worldNodeIndices.size());
                 ++localIndex)
            {
                if (worldNodeIndices[static_cast<std::size_t>(localIndex)] == worldNodeIndex)
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
    }

    void TetMeshComponent::SetGeometry(
        const Vec3* positions,
        int nodeCount,
        const int* tetLocalNodeIndices,
        int tetCount)
    {
        restNodePositions.clear();
        currentNodePositions.clear();
        nodeIndices.clear();
        tets.clear();

        if (positions != nullptr && nodeCount > 0)
        {
            restNodePositions.assign(positions, positions + nodeCount);
            currentNodePositions = restNodePositions;
            nodeIndices.reserve(static_cast<std::size_t>(nodeCount));
            for (int i = 0; i < nodeCount; ++i)
            {
                nodeIndices.push_back(i);
            }
        }

        if (tetLocalNodeIndices != nullptr && tetCount > 0)
        {
            tets.reserve(static_cast<std::size_t>(tetCount));
            for (int tetIndex = 0; tetIndex < tetCount; ++tetIndex)
            {
                const int local0 = tetLocalNodeIndices[tetIndex * 4 + 0];
                const int local1 = tetLocalNodeIndices[tetIndex * 4 + 1];
                const int local2 = tetLocalNodeIndices[tetIndex * 4 + 2];
                const int local3 = tetLocalNodeIndices[tetIndex * 4 + 3];
                const int localNodeCount = static_cast<int>(restNodePositions.size());

                if (local0 < 0 || local0 >= localNodeCount ||
                    local1 < 0 || local1 >= localNodeCount ||
                    local2 < 0 || local2 >= localNodeCount ||
                    local3 < 0 || local3 >= localNodeCount)
                {
                    continue;
                }

                tets.push_back(MakeTet(local0, local1, local2, local3));
            }
        }
    }

    int TetMeshComponent::GetNodeCount() const
    {
        return static_cast<int>(restNodePositions.size());
    }

    int TetMeshComponent::GetTetCount() const
    {
        return static_cast<int>(tets.size());
    }

    int TetMeshComponent::GetTetNodeIndex(int tetIndex, int cornerIndex) const
    {
        if (tetIndex < 0 || tetIndex >= static_cast<int>(tets.size()))
        {
            return -1;
        }

        const Tet& tet = tets[static_cast<std::size_t>(tetIndex)];
        switch (cornerIndex)
        {
        case 0:
            return tet.node0;
        case 1:
            return tet.node1;
        case 2:
            return tet.node2;
        case 3:
            return tet.node3;
        default:
            return -1;
        }
    }

    const Vec3& TetMeshComponent::GetLocalRestPosition(int localNodeIndex) const
    {
        if (localNodeIndex < 0 ||
            localNodeIndex >= static_cast<int>(restNodePositions.size()))
        {
            return ZeroVec3;
        }

        return restNodePositions[static_cast<std::size_t>(localNodeIndex)];
    }

    const Vec3& TetMeshComponent::GetLocalCurrentPosition(int localNodeIndex) const
    {
        if (localNodeIndex < 0 ||
            localNodeIndex >= static_cast<int>(currentNodePositions.size()))
        {
            return ZeroVec3;
        }

        return currentNodePositions[static_cast<std::size_t>(localNodeIndex)];
    }

    int TetMeshComponent::GetWorldNodeIndex(int localNodeIndex) const
    {
        (void)localNodeIndex;
        return -1;
    }

    void TetMeshComponent::SetLocalCurrentPosition(
        int localNodeIndex,
        const Vec3& position)
    {
        if (localNodeIndex < 0 ||
            localNodeIndex >= static_cast<int>(currentNodePositions.size()))
        {
            return;
        }

        currentNodePositions[static_cast<std::size_t>(localNodeIndex)] = position;
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
        topologyDirty = true;
    }

    void TetMeshComponent::DeactivateTet(int tetIndex)
    {
        SetTetActive(tetIndex, false);
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

    void TetMeshComponent::PostUpdate(World& world, float dt)
    {
        (void)dt;

        if (!topologyDirty)
        {
            return;
        }

        PhysicsEvent event;
        event.type = PhysicsEventType::TetMeshTopologyChanged;
        event.world = &world;
        event.sender = this;
        world.EmitEvent(event);

        topologyDirty = false;
    }

    std::unique_ptr<TetMeshPhysicsComponent> TetMeshPhysicsComponent::CreateFromGlobalNodes(
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

    std::unique_ptr<TetMeshPhysicsComponent> TetMeshPhysicsComponent::CreateFromGlobalNodes(
        World& world,
        const int* globalNodeIndices,
        int nodeCount,
        const int* tetGlobalNodeIndices,
        int tetCount,
        const TetMeshComponentDesc& desc)
    {
        auto component = std::make_unique<TetMeshPhysicsComponent>();
        component->material = desc.material;
        component->selectedFemModel = desc.femModel;

        if (globalNodeIndices != nullptr && nodeCount > 0)
        {
            component->worldNodeIndices.assign(
                globalNodeIndices,
                globalNodeIndices + nodeCount);
            component->restNodePositions.reserve(static_cast<std::size_t>(nodeCount));
            component->currentNodePositions.reserve(static_cast<std::size_t>(nodeCount));
            component->nodeIndices.reserve(static_cast<std::size_t>(nodeCount));
            for (int i = 0; i < nodeCount; ++i)
            {
                const int worldNodeIndex = globalNodeIndices[i];
                if (worldNodeIndex < 0 ||
                    worldNodeIndex >= static_cast<int>(world.GetNodes().size()))
                {
                    continue;
                }

                const Vec3& position =
                    world.GetNode(worldNodeIndex).restPosition;
                component->restNodePositions.push_back(position);
                component->currentNodePositions.push_back(
                    world.GetNode(worldNodeIndex).position);
                component->nodeIndices.push_back(
                    static_cast<int>(component->nodeIndices.size()));
            }
        }

        if (tetGlobalNodeIndices != nullptr && tetCount > 0)
        {
            component->tets.reserve(static_cast<std::size_t>(tetCount));
            for (int i = 0; i < tetCount; ++i)
            {
                const int local0 = FindLocalNodeIndex(
                    component->worldNodeIndices,
                    tetGlobalNodeIndices[i * 4 + 0]);
                const int local1 = FindLocalNodeIndex(
                    component->worldNodeIndices,
                    tetGlobalNodeIndices[i * 4 + 1]);
                const int local2 = FindLocalNodeIndex(
                    component->worldNodeIndices,
                    tetGlobalNodeIndices[i * 4 + 2]);
                const int local3 = FindLocalNodeIndex(
                    component->worldNodeIndices,
                    tetGlobalNodeIndices[i * 4 + 3]);

                if (local0 < 0 || local1 < 0 || local2 < 0 || local3 < 0)
                {
                    continue;
                }

                component->tets.push_back(MakeTet(local0, local1, local2, local3));
            }
        }

        component->SyncWorldTetsFromLocalTets(world);
        component->RebuildTetFemCache();
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

    std::unique_ptr<TetMeshPhysicsComponent> TetMeshPhysicsComponent::CreateFromPositions(
        World& world,
        const Vec3* positions,
        const int* fixedNodeFlags,
        int nodeCount,
        const int* tetLocalNodeIndices,
        int tetCount,
        const TetMeshComponentDesc& desc)
    {
        auto component = std::make_unique<TetMeshPhysicsComponent>();
        component->material = desc.material;
        component->selectedFemModel = desc.femModel;

        if (positions != nullptr && nodeCount > 0)
        {
            component->worldNodeIndices.reserve(static_cast<std::size_t>(nodeCount));
            for (int i = 0; i < nodeCount; ++i)
            {
                const int nodeIndex = world.AddNode(positions[i]);
                if (fixedNodeFlags != nullptr && fixedNodeFlags[i] != 0)
                {
                    world.SetNodeFixed(nodeIndex, true);
                }
                component->worldNodeIndices.push_back(nodeIndex);
            }
        }

        component->SetGeometry(positions, nodeCount, tetLocalNodeIndices, tetCount);
        component->SyncWorldTetsFromLocalTets(world);
        component->RebuildTetFemCache();
        return component;
    }

    void TetMeshPhysicsComponent::SetMaterial(const Material& value)
    {
        material = value;
        for (Tet& tet : worldTets)
        {
            ApplyMaterialToTet(tet, material);
        }
        RebuildTetFemCache();
    }

    void TetMeshPhysicsComponent::SyncCurrentPositionsFromWorld(const World& world)
    {
        currentNodePositions.resize(worldNodeIndices.size());
        for (int localIndex = 0;
             localIndex < static_cast<int>(worldNodeIndices.size());
             ++localIndex)
        {
            const int worldNodeIndex = worldNodeIndices[static_cast<std::size_t>(localIndex)];
            if (worldNodeIndex < 0 ||
                worldNodeIndex >= static_cast<int>(world.GetNodes().size()))
            {
                continue;
            }

            currentNodePositions[static_cast<std::size_t>(localIndex)] =
                world.GetNode(worldNodeIndex).position;
        }
    }

    void TetMeshPhysicsComponent::SyncWorldTetsFromLocalTets(const World& world)
    {
        (void)world;
        worldTets.clear();
        worldTets.reserve(tets.size());
        for (const Tet& localTet : tets)
        {
            const int localNodes[4] = {
                localTet.node0,
                localTet.node1,
                localTet.node2,
                localTet.node3};
            bool valid = true;
            int worldNodes[4] = {};
            for (int corner = 0; corner < 4; ++corner)
            {
                const int localNode = localNodes[corner];
                if (localNode < 0 ||
                    localNode >= static_cast<int>(worldNodeIndices.size()))
                {
                    valid = false;
                    break;
                }

                worldNodes[corner] =
                    worldNodeIndices[static_cast<std::size_t>(localNode)];
            }

            if (!valid)
            {
                continue;
            }

            Tet worldTet = MakeTet(
                worldNodes[0],
                worldNodes[1],
                worldNodes[2],
                worldNodes[3]);
            worldTet.active = localTet.active;
            ApplyMaterialToTet(worldTet, material);
            FEMModel::InitializeTetRestData(worldTet, world.GetNodes());
            worldTets.push_back(worldTet);
        }
    }

    void TetMeshPhysicsComponent::RebuildTetFemCache()
    {
        tetFemCache.clear();
        tetFemCache.reserve(worldTets.size());
        for (const Tet& tet : worldTets)
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
            FEMModel::BuildSparsePatternFromTetConnectivity(worldTets));
        femSparsePatternDirty = false;
    }

    int TetMeshPhysicsComponent::GetWorldNodeIndex(int localNodeIndex) const
    {
        if (localNodeIndex < 0 ||
            localNodeIndex >= static_cast<int>(worldNodeIndices.size()))
        {
            return -1;
        }

        return worldNodeIndices[static_cast<std::size_t>(localNodeIndex)];
    }

    void TetMeshPhysicsComponent::SetLocalCurrentPosition(
        int localNodeIndex,
        const Vec3& position)
    {
        TetMeshComponent::SetLocalCurrentPosition(localNodeIndex, position);
    }

    void TetMeshPhysicsComponent::UpdateSystem(
        World& world,
        SolverData& solverData,
        float dt)
    {
        (void)dt;
        if (!physicsEnabled)
        {
            SyncCurrentPositionsFromWorld(world);
            return;
        }

        SyncCurrentPositionsFromWorld(world);
        SyncWorldTetsFromLocalTets(world);
        FEMModel::AssembleLumpedMass(*this, world, solverData);
        EnsureFemSparsePattern(static_cast<int>(world.GetNodes().size()));
        if (tetFemCache.size() != worldTets.size())
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
