#include "PhysiK/Components/TetMeshMapperComponent.h"

#include <cstddef>

#include "PhysiK/Components/TetMeshComponent.h"
#include "PhysiK/Components/VisualMeshComponent.h"
#include "PhysiK/Core/World/World.h"
#include "PhysiK/PhysicsData/Node.h"
#include "PhysiK/PhysicsData/Tet.h"

namespace PhysiK
{
    TetMeshMapperComponent::TetMeshMapperComponent(
        ComponentHandle sourceTetMeshHandle,
        ComponentHandle destinationTetMeshHandle)
    {
        this->sourceTetMeshHandle = sourceTetMeshHandle;
        this->destinationTetMeshHandle = destinationTetMeshHandle;
    }

    void TetMeshMapperComponent::BuildTetMeshMapping(const World& world)
    {
        embeddedDestinationVertices.clear();

        const Component* sourceComponent = world.GetComponent(sourceTetMeshHandle);
        const TetMeshComponent* sourceTetMesh =
            dynamic_cast<const TetMeshComponent*>(sourceComponent);
        const Component* destinationComponent =
            world.GetComponent(destinationTetMeshHandle);
        const TetMeshComponent* destinationTetMesh =
            dynamic_cast<const TetMeshComponent*>(destinationComponent);
        if (sourceTetMesh == nullptr || destinationTetMesh == nullptr)
        {
            return;
        }

        const std::vector<Node>& nodes = world.GetNodes();
        embeddedDestinationVertices.resize(destinationTetMesh->nodeIndices.size());
        constexpr float insideEpsilon = -0.0001f;

        for (std::size_t destinationVertex = 0;
             destinationVertex < destinationTetMesh->nodeIndices.size();
             ++destinationVertex)
        {
            TetMeshMappedVertex& mappedVertex =
                embeddedDestinationVertices[destinationVertex];
            mappedVertex = TetMeshMappedVertex{};

            const int destinationNode =
                destinationTetMesh->nodeIndices[destinationVertex];
            if (destinationNode < 0 ||
                destinationNode >= static_cast<int>(nodes.size()))
            {
                continue;
            }

            const Vec3& destinationRestPosition =
                nodes[static_cast<std::size_t>(destinationNode)].restPosition;

            for (std::size_t sourceTetIndex = 0;
                 sourceTetIndex < sourceTetMesh->tets.size();
                 ++sourceTetIndex)
            {
                const Tet& sourceTet = sourceTetMesh->tets[sourceTetIndex];
                if (!sourceTet.active)
                {
                    continue;
                }

                const int sourceNodes[4] = {
                    sourceTet.node0,
                    sourceTet.node1,
                    sourceTet.node2,
                    sourceTet.node3};
                bool sourceNodesAreValid = true;
                for (int sourceNode : sourceNodes)
                {
                    if (sourceNode < 0 ||
                        sourceNode >= static_cast<int>(nodes.size()))
                    {
                        sourceNodesAreValid = false;
                        break;
                    }
                }

                if (!sourceNodesAreValid)
                {
                    continue;
                }

                Vec4 barycentric;
                if (!ComputeTetBarycentric(
                        destinationRestPosition,
                        nodes[static_cast<std::size_t>(sourceTet.node0)].restPosition,
                        nodes[static_cast<std::size_t>(sourceTet.node1)].restPosition,
                        nodes[static_cast<std::size_t>(sourceTet.node2)].restPosition,
                        nodes[static_cast<std::size_t>(sourceTet.node3)].restPosition,
                        barycentric))
                {
                    continue;
                }

                if (barycentric.x >= insideEpsilon &&
                    barycentric.y >= insideEpsilon &&
                    barycentric.z >= insideEpsilon &&
                    barycentric.w >= insideEpsilon)
                {
                    mappedVertex.sourceTetIndex =
                        static_cast<int>(sourceTetIndex);
                    mappedVertex.barycentric = barycentric;
                    mappedVertex.valid = true;
                    break;
                }
            }
        }
    }

    void TetMeshMapperComponent::UpdateDestinationNodes(World& world)
    {
        const Component* sourceComponent = world.GetComponent(sourceTetMeshHandle);
        const TetMeshComponent* sourceTetMesh =
            dynamic_cast<const TetMeshComponent*>(sourceComponent);
        Component* destinationComponent =
            world.GetComponent(destinationTetMeshHandle);
        TetMeshComponent* destinationTetMesh =
            dynamic_cast<TetMeshComponent*>(destinationComponent);
        if (sourceTetMesh == nullptr || destinationTetMesh == nullptr)
        {
            return;
        }

        if (embeddedDestinationVertices.size() !=
            destinationTetMesh->nodeIndices.size())
        {
            return;
        }

        const std::vector<Node>& nodes = world.GetNodes();
        for (std::size_t destinationVertex = 0;
             destinationVertex < embeddedDestinationVertices.size();
             ++destinationVertex)
        {
            const TetMeshMappedVertex& mappedVertex =
                embeddedDestinationVertices[destinationVertex];
            if (!mappedVertex.valid ||
                mappedVertex.sourceTetIndex < 0 ||
                mappedVertex.sourceTetIndex >=
                    static_cast<int>(sourceTetMesh->tets.size()))
            {
                continue;
            }

            const Tet& sourceTet = sourceTetMesh->tets[
                static_cast<std::size_t>(mappedVertex.sourceTetIndex)];
            if (!sourceTet.active)
            {
                continue;
            }

            const int sourceNodes[4] = {
                sourceTet.node0,
                sourceTet.node1,
                sourceTet.node2,
                sourceTet.node3};
            bool sourceNodesAreValid = true;
            for (int sourceNode : sourceNodes)
            {
                if (sourceNode < 0 || sourceNode >= static_cast<int>(nodes.size()))
                {
                    sourceNodesAreValid = false;
                    break;
                }
            }

            const int destinationNode =
                destinationTetMesh->nodeIndices[destinationVertex];
            if (!sourceNodesAreValid ||
                destinationNode < 0 ||
                destinationNode >= static_cast<int>(nodes.size()))
            {
                continue;
            }

            const Vec4& weights = mappedVertex.barycentric;
            Node& destination = world.GetNode(destinationNode);
            destination.position =
                nodes[static_cast<std::size_t>(sourceTet.node0)].position * weights.x +
                nodes[static_cast<std::size_t>(sourceTet.node1)].position * weights.y +
                nodes[static_cast<std::size_t>(sourceTet.node2)].position * weights.z +
                nodes[static_cast<std::size_t>(sourceTet.node3)].position * weights.w;
            destination.velocity = Vec3{};
        }
    }

    void TetMeshMapperComponent::PostUpdate(World& world, float dt)
    {
        (void)dt;
        UpdateDestinationNodes(world);
    }
}
