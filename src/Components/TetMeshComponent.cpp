#include "PhysiK/Components/TetMeshComponent.h"

#include <algorithm>
#include <cmath>

#include "PhysiK/Core/Solvers/SolverData.h"
#include "PhysiK/Core/World/World.h"

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
        }

        void AddLumpedMassToNode(World& world, int nodeIndex, float mass)
        {
            if (nodeIndex < 0 || !std::isfinite(mass) || mass <= 0.0f)
            {
                return;
            }

            Node& node = world.GetNode(nodeIndex);
            if (node.inverseMass <= 0.0f)
            {
                return;
            }

            node.femMass += mass;
            if (node.femMass > 0.0f && std::isfinite(node.femMass))
            {
                node.inverseMass = 1.0f / node.femMass;
            }
        }

        void ApplyLumpedMass(TetMeshComponent& component, World& world)
        {
            const float density = std::max(0.0f, component.material.density);
            if (!std::isfinite(density) || density <= 0.0f)
            {
                return;
            }

            for (const Tet& tet : component.tets)
            {
                if (!std::isfinite(tet.restVolume) || tet.restVolume <= 0.0f)
                {
                    continue;
                }

                const float nodalMass = density * tet.restVolume * 0.25f;
                AddLumpedMassToNode(world, tet.node0, nodalMass);
                AddLumpedMassToNode(world, tet.node1, nodalMass);
                AddLumpedMassToNode(world, tet.node2, nodalMass);
                AddLumpedMassToNode(world, tet.node3, nodalMass);
            }
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
        auto component = std::make_unique<TetMeshComponent>();
        component->material = material;

        if (globalNodeIndices != nullptr && nodeCount > 0)
        {
            component->nodeIndices.assign(globalNodeIndices, globalNodeIndices + nodeCount);
        }

        if (tetGlobalNodeIndices != nullptr && tetCount > 0)
        {
            component->tets.reserve(static_cast<std::size_t>(tetCount));
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

        ApplyLumpedMass(*component, world);

        return component;
    }

    std::unique_ptr<TetMeshComponent> TetMeshComponent::CreateFromPositions(
        World& world,
        const Vec3* positions,
        const float* inverseMasses,
        int nodeCount,
        const int* tetLocalNodeIndices,
        int tetCount,
        const Material& material)
    {
        auto component = std::make_unique<TetMeshComponent>();
        component->material = material;

        if (positions != nullptr && nodeCount > 0)
        {
            component->nodeIndices.reserve(static_cast<std::size_t>(nodeCount));
            for (int i = 0; i < nodeCount; ++i)
            {
                const int nodeIndex = world.AddNode(positions[i]);
                if (inverseMasses != nullptr && inverseMasses[i] <= 0.0f)
                {
                    world.SetNodeFixed(nodeIndex, true);
                }
                component->nodeIndices.push_back(nodeIndex);
            }
        }

        if (tetLocalNodeIndices != nullptr && tetCount > 0)
        {
            component->tets.reserve(static_cast<std::size_t>(tetCount));
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

        ApplyLumpedMass(*component, world);

        return component;
    }

    void TetMeshComponent::UpdateSystem(
        World& world,
        SolverData& solverData,
        float dt)
    {
        femModel.UpdateSystem(world, *this, solverData, dt);
    }
}
