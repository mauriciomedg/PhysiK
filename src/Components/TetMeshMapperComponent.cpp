#include "PhysiK/Components/TetMeshMapperComponent.h"

#include <cstddef>

#include "PhysiK/Components/TetMeshComponent.h"
#include "PhysiK/Components/VisualMeshComponent.h"
#include "PhysiK/Core/World/World.h"
#include "PhysiK/PhysicsData/Tet.h"

namespace PhysiK
{
    TetMeshMapperComponent::TetMeshMapperComponent()
    {
        listenedEvents.push_back(PhysicsEventType::TetMeshTopologyChanged);
    }

    TetMeshMapperComponent::TetMeshMapperComponent(
        ComponentHandle sourceTetMeshHandle,
        ComponentHandle destinationTetMeshHandle)
        : sourceTetMeshHandle(sourceTetMeshHandle)
        , destinationTetMeshHandle(destinationTetMeshHandle)
        , mappingDirty(true)
    {
        listenedEvents.push_back(PhysicsEventType::TetMeshTopologyChanged);
    }

    bool TetMeshMapperComponent::BuildTetMeshMapping(World& world)
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
            return false;
        }

        if (sourceTetMesh->GetTetCount() <= 0 ||
            destinationTetMesh->GetNodeCount() <= 0)
        {
            return false;
        }

        embeddedDestinationVertices.resize(
            static_cast<std::size_t>(destinationTetMesh->GetNodeCount()));
        constexpr float insideEpsilon = -0.0001f;

        for (std::size_t destinationVertex = 0;
             destinationVertex < embeddedDestinationVertices.size();
             ++destinationVertex)
        {
            TetMeshMappedVertex& mappedVertex =
                embeddedDestinationVertices[destinationVertex];
            mappedVertex = TetMeshMappedVertex{};

            const int destinationNode = static_cast<int>(destinationVertex);
            if (destinationNode < 0 ||
                destinationNode >= destinationTetMesh->GetNodeCount())
            {
                continue;
            }

            const Vec3& destinationRestPosition =
                destinationTetMesh->GetLocalRestPosition(destinationNode);

            for (int sourceTetIndex = 0;
                 sourceTetIndex < sourceTetMesh->GetTetCount();
                 ++sourceTetIndex)
            {
                if (!sourceTetMesh->IsTetActive(sourceTetIndex))
                {
                    continue;
                }

                const int sourceNodes[4] = {
                    sourceTetMesh->GetTetNodeIndex(sourceTetIndex, 0),
                    sourceTetMesh->GetTetNodeIndex(sourceTetIndex, 1),
                    sourceTetMesh->GetTetNodeIndex(sourceTetIndex, 2),
                    sourceTetMesh->GetTetNodeIndex(sourceTetIndex, 3)};
                bool sourceNodesAreValid = true;
                for (int sourceNode : sourceNodes)
                {
                    if (sourceNode < 0 ||
                        sourceNode >= sourceTetMesh->GetNodeCount())
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
                        sourceTetMesh->GetLocalRestPosition(sourceNodes[0]),
                        sourceTetMesh->GetLocalRestPosition(sourceNodes[1]),
                        sourceTetMesh->GetLocalRestPosition(sourceNodes[2]),
                        sourceTetMesh->GetLocalRestPosition(sourceNodes[3]),
                        barycentric))
                {
                    continue;
                }

                if (barycentric.x >= insideEpsilon &&
                    barycentric.y >= insideEpsilon &&
                    barycentric.z >= insideEpsilon &&
                    barycentric.w >= insideEpsilon)
                {
                    mappedVertex.sourceTetIndex = sourceTetIndex;
                    mappedVertex.barycentric = barycentric;
                    mappedVertex.valid = true;
                    break;
                }
            }
        }

        return true;
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
            static_cast<std::size_t>(destinationTetMesh->GetNodeCount()))
        {
            return;
        }

        for (std::size_t destinationVertex = 0;
             destinationVertex < embeddedDestinationVertices.size();
             ++destinationVertex)
        {
            const TetMeshMappedVertex& mappedVertex =
                embeddedDestinationVertices[destinationVertex];
            if (!mappedVertex.valid ||
                mappedVertex.sourceTetIndex < 0 ||
                mappedVertex.sourceTetIndex >=
                    sourceTetMesh->GetTetCount())
            {
                continue;
            }

            if (!sourceTetMesh->IsTetActive(mappedVertex.sourceTetIndex))
            {
                continue;
            }

            const int sourceNodes[4] = {
                sourceTetMesh->GetTetNodeIndex(mappedVertex.sourceTetIndex, 0),
                sourceTetMesh->GetTetNodeIndex(mappedVertex.sourceTetIndex, 1),
                sourceTetMesh->GetTetNodeIndex(mappedVertex.sourceTetIndex, 2),
                sourceTetMesh->GetTetNodeIndex(mappedVertex.sourceTetIndex, 3)};
            bool sourceNodesAreValid = true;
            for (int sourceNode : sourceNodes)
            {
                if (sourceNode < 0 || sourceNode >= sourceTetMesh->GetNodeCount())
                {
                    sourceNodesAreValid = false;
                    break;
                }
            }

            const int destinationNode =
                static_cast<int>(destinationVertex);
            if (!sourceNodesAreValid ||
                destinationNode < 0 ||
                destinationNode >= destinationTetMesh->GetNodeCount())
            {
                continue;
            }

            const Vec4& weights = mappedVertex.barycentric;
            const Vec3 mappedPosition =
                sourceTetMesh->GetLocalCurrentPosition(sourceNodes[0]) * weights.x +
                sourceTetMesh->GetLocalCurrentPosition(sourceNodes[1]) * weights.y +
                sourceTetMesh->GetLocalCurrentPosition(sourceNodes[2]) * weights.z +
                sourceTetMesh->GetLocalCurrentPosition(sourceNodes[3]) * weights.w;
            destinationTetMesh->SetLocalCurrentPosition(destinationNode, mappedPosition);
        }
    }

    void TetMeshMapperComponent::MarkMappingDirty()
    {
        mappingDirty = true;
    }

    bool TetMeshMapperComponent::IsMappingDirty() const
    {
        return mappingDirty;
    }

    void TetMeshMapperComponent::OnPhysicsEvent(const PhysicsEvent& event)
    {
        if (event.type == PhysicsEventType::TetMeshTopologyChanged)
        {
            MarkMappingDirty();
        }
    }

    void TetMeshMapperComponent::PostUpdate(World& world, float dt)
    {
        (void)dt;
        if (mappingDirty)
        {
            if (!BuildTetMeshMapping(world))
            {
                return;
            }

            mappingDirty = false;
        }

        UpdateDestinationNodes(world);
    }
}
