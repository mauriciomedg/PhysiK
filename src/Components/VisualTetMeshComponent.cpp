#include "PhysiK/Components/VisualTetMeshComponent.h"

#include <cstddef>

#include "PhysiK/Components/TetMeshComponent.h"
#include "PhysiK/Components/VisualMeshComponent.h"
#include "PhysiK/Core/World/World.h"
#include "PhysiK/PhysicsData/Node.h"
#include "PhysiK/PhysicsData/Tet.h"

namespace PhysiK
{
    VisualTetMeshComponent::VisualTetMeshComponent(
        ComponentHandle hostMechanicalTetMeshHandle)
    {
        this->hostMechanicalTetMeshHandle = hostMechanicalTetMeshHandle;
    }

    void VisualTetMeshComponent::SetVisualTetMeshData(
        const Vec3* vertices,
        int vertexCount,
        const int* tetIndices,
        int tetIndexCount)
    {
        restVisualVertices.clear();
        deformedVisualVertices.clear();
        visualTetIndices.clear();
        embeddedVertices.clear();
        activeVisualTets.clear();

        if (vertices != nullptr && vertexCount > 0)
        {
            restVisualVertices.assign(vertices, vertices + vertexCount);
            deformedVisualVertices = restVisualVertices;
            embeddedVertices.resize(static_cast<std::size_t>(vertexCount));
        }

        if (tetIndices != nullptr && tetIndexCount > 0)
        {
            visualTetIndices.assign(tetIndices, tetIndices + tetIndexCount);
            activeVisualTets.assign(
                static_cast<std::size_t>(tetIndexCount / 4),
                true);
        }
    }

    void VisualTetMeshComponent::BuildEmbedding(const World& world)
    {
        embeddedVertices.clear();
        embeddedVertices.resize(restVisualVertices.size());

        const Component* hostComponent =
            world.GetComponent(hostMechanicalTetMeshHandle);
        const TetMeshComponent* hostMechanicalTetMesh =
            dynamic_cast<const TetMeshComponent*>(hostComponent);
        if (hostMechanicalTetMesh == nullptr)
        {
            return;
        }

        const std::vector<Node>& nodes = world.GetNodes();
        constexpr float insideEpsilon = -0.0001f;

        for (std::size_t vertexIndex = 0;
             vertexIndex < restVisualVertices.size();
             ++vertexIndex)
        {
            VisualTetEmbeddedVertex& embeddedVertex =
                embeddedVertices[vertexIndex];
            embeddedVertex = VisualTetEmbeddedVertex{};

            for (std::size_t tetIndex = 0;
                 tetIndex < hostMechanicalTetMesh->tets.size();
                 ++tetIndex)
            {
                const Tet& tet = hostMechanicalTetMesh->tets[tetIndex];
                if (!tet.active)
                {
                    continue;
                }

                const int nodeIndices[4] = {tet.node0, tet.node1, tet.node2, tet.node3};
                bool nodesAreValid = true;
                for (int nodeIndex : nodeIndices)
                {
                    if (nodeIndex < 0 ||
                        nodeIndex >= static_cast<int>(nodes.size()))
                    {
                        nodesAreValid = false;
                        break;
                    }
                }

                if (!nodesAreValid)
                {
                    continue;
                }

                Vec4 barycentric;
                if (!ComputeTetBarycentric(
                        restVisualVertices[vertexIndex],
                        nodes[static_cast<std::size_t>(tet.node0)].restPosition,
                        nodes[static_cast<std::size_t>(tet.node1)].restPosition,
                        nodes[static_cast<std::size_t>(tet.node2)].restPosition,
                        nodes[static_cast<std::size_t>(tet.node3)].restPosition,
                        barycentric))
                {
                    continue;
                }

                if (barycentric.x >= insideEpsilon &&
                    barycentric.y >= insideEpsilon &&
                    barycentric.z >= insideEpsilon &&
                    barycentric.w >= insideEpsilon)
                {
                    embeddedVertex.mechanicalTetIndex = static_cast<int>(tetIndex);
                    embeddedVertex.barycentric = barycentric;
                    embeddedVertex.valid = true;
                    break;
                }
            }
        }
    }

    void VisualTetMeshComponent::UpdateDeformedVertices(const World& world)
    {
        if (deformedVisualVertices.size() != restVisualVertices.size())
        {
            deformedVisualVertices = restVisualVertices;
        }

        if (embeddedVertices.size() != restVisualVertices.size())
        {
            return;
        }

        const Component* hostComponent =
            world.GetComponent(hostMechanicalTetMeshHandle);
        const TetMeshComponent* hostMechanicalTetMesh =
            dynamic_cast<const TetMeshComponent*>(hostComponent);
        if (hostMechanicalTetMesh == nullptr)
        {
            return;
        }

        const std::vector<Node>& nodes = world.GetNodes();

        for (std::size_t vertexIndex = 0;
             vertexIndex < embeddedVertices.size();
             ++vertexIndex)
        {
            const VisualTetEmbeddedVertex& embeddedVertex =
                embeddedVertices[vertexIndex];
            if (!embeddedVertex.valid ||
                embeddedVertex.mechanicalTetIndex < 0 ||
                embeddedVertex.mechanicalTetIndex >=
                    static_cast<int>(hostMechanicalTetMesh->tets.size()))
            {
                continue;
            }

            const Tet& tet = hostMechanicalTetMesh->tets[
                static_cast<std::size_t>(embeddedVertex.mechanicalTetIndex)];
            if (!tet.active)
            {
                continue;
            }

            const int nodeIndices[4] = {tet.node0, tet.node1, tet.node2, tet.node3};
            bool nodesAreValid = true;
            for (int nodeIndex : nodeIndices)
            {
                if (nodeIndex < 0 || nodeIndex >= static_cast<int>(nodes.size()))
                {
                    nodesAreValid = false;
                    break;
                }
            }

            if (!nodesAreValid)
            {
                continue;
            }

            const Vec4& weights = embeddedVertex.barycentric;
            deformedVisualVertices[vertexIndex] =
                nodes[static_cast<std::size_t>(tet.node0)].position * weights.x +
                nodes[static_cast<std::size_t>(tet.node1)].position * weights.y +
                nodes[static_cast<std::size_t>(tet.node2)].position * weights.z +
                nodes[static_cast<std::size_t>(tet.node3)].position * weights.w;
        }
    }

    const std::vector<Vec3>& VisualTetMeshComponent::GetDeformedVertices() const
    {
        return deformedVisualVertices;
    }

    const std::vector<int>& VisualTetMeshComponent::GetVisualTetIndices() const
    {
        return visualTetIndices;
    }

    void VisualTetMeshComponent::PostUpdate(World& world, float dt)
    {
        (void)dt;
        UpdateDeformedVertices(world);
    }
}
