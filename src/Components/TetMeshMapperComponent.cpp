#include "PhysiK/Components/TetMeshMapperComponent.h"

#include <cstddef>

#include "PhysiK/Components/TetMeshComponent.h"
#include "PhysiK/Components/VisualMeshComponent.h"
#include "PhysiK/Core/World/World.h"
#include "PhysiK/PhysicsData/Tet.h"

namespace PhysiK
{
    namespace
    {
        bool IsLocalTetNodeValid(const TetMeshComponent& mesh, int localNodeIndex)
        {
            return localNodeIndex >= 0 && localNodeIndex < mesh.GetNodeCount();
        }

        bool TryGetTetNodes(
            const TetMeshComponent& mesh,
            int tetIndex,
            int (&outNodes)[4])
        {
            if (tetIndex < 0 || tetIndex >= mesh.GetTetCount())
            {
                return false;
            }

            outNodes[0] = mesh.GetTetNodeIndex(tetIndex, 0);
            outNodes[1] = mesh.GetTetNodeIndex(tetIndex, 1);
            outNodes[2] = mesh.GetTetNodeIndex(tetIndex, 2);
            outNodes[3] = mesh.GetTetNodeIndex(tetIndex, 3);

            for (int nodeIndex : outNodes)
            {
                if (!IsLocalTetNodeValid(mesh, nodeIndex))
                {
                    return false;
                }
            }

            return true;
        }

        int FindContainingActiveSourceTet(
            const TetMeshComponent& sourceTetMesh,
            const Vec3& point)
        {
            constexpr float insideEpsilon = -0.0001f;

            for (int sourceTetIndex = 0;
                 sourceTetIndex < sourceTetMesh.GetTetCount();
                 ++sourceTetIndex)
            {
                if (!sourceTetMesh.IsTetActive(sourceTetIndex))
                {
                    continue;
                }

                int sourceNodes[4];
                if (!TryGetTetNodes(sourceTetMesh, sourceTetIndex, sourceNodes))
                {
                    continue;
                }

                Vec4 barycentric;
                if (!ComputeTetBarycentric(
                        point,
                        sourceTetMesh.GetLocalRestPosition(sourceNodes[0]),
                        sourceTetMesh.GetLocalRestPosition(sourceNodes[1]),
                        sourceTetMesh.GetLocalRestPosition(sourceNodes[2]),
                        sourceTetMesh.GetLocalRestPosition(sourceNodes[3]),
                        barycentric))
                {
                    continue;
                }

                if (barycentric.x >= insideEpsilon &&
                    barycentric.y >= insideEpsilon &&
                    barycentric.z >= insideEpsilon &&
                    barycentric.w >= insideEpsilon)
                {
                    return sourceTetIndex;
                }
            }

            return -1;
        }
    }

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
        , activeStateDirty(true)
    {
        listenedEvents.push_back(PhysicsEventType::TetMeshTopologyChanged);
    }

    bool TetMeshMapperComponent::BuildTetMeshMapping(World& world)
    {
        embeddedDestinationVertices.clear();
        embeddedDestinationTets.clear();

        const Component* sourceComponent = world.GetComponent(sourceTetMeshHandle);
        const TetMeshComponent* sourceTetMesh =
            dynamic_cast<const TetMeshComponent*>(sourceComponent);
        Component* destinationComponent =
            world.GetComponent(destinationTetMeshHandle);
        TetMeshComponent* destinationTetMesh =
            dynamic_cast<TetMeshComponent*>(destinationComponent);
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
        embeddedDestinationTets.resize(
            static_cast<std::size_t>(destinationTetMesh->GetTetCount()));
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

                int sourceNodes[4];
                if (!TryGetTetNodes(*sourceTetMesh, sourceTetIndex, sourceNodes))
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

        for (int destinationTetIndex = 0;
             destinationTetIndex < destinationTetMesh->GetTetCount();
             ++destinationTetIndex)
        {
            TetMeshMappedTet& mappedTet =
                embeddedDestinationTets[static_cast<std::size_t>(destinationTetIndex)];
            mappedTet = TetMeshMappedTet{};

            int destinationNodes[4];
            if (!TryGetTetNodes(
                    *destinationTetMesh,
                    destinationTetIndex,
                    destinationNodes))
            {
                continue;
            }

            const Vec3 centroid =
                (destinationTetMesh->GetLocalRestPosition(destinationNodes[0]) +
                 destinationTetMesh->GetLocalRestPosition(destinationNodes[1]) +
                 destinationTetMesh->GetLocalRestPosition(destinationNodes[2]) +
                 destinationTetMesh->GetLocalRestPosition(destinationNodes[3])) *
                0.25f;

            const int sourceTetIndex =
                FindContainingActiveSourceTet(*sourceTetMesh, centroid);
            if (sourceTetIndex >= 0)
            {
                mappedTet.sourceTetIndex = sourceTetIndex;
                mappedTet.valid = true;
            }
        }

        return true;
    }

    bool TetMeshMapperComponent::RefreshDestinationActiveStates(World& world)
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
            return false;
        }

        if (embeddedDestinationTets.size() !=
                static_cast<std::size_t>(destinationTetMesh->GetTetCount()) ||
            embeddedDestinationVertices.size() !=
                static_cast<std::size_t>(destinationTetMesh->GetNodeCount()))
        {
            return false;
        }

        for (int destinationTetIndex = 0;
             destinationTetIndex < destinationTetMesh->GetTetCount();
             ++destinationTetIndex)
        {
            const TetMeshMappedTet& mappedTet =
                embeddedDestinationTets[static_cast<std::size_t>(destinationTetIndex)];
            bool shouldBeActive =
                mappedTet.valid &&
                mappedTet.sourceTetIndex >= 0 &&
                mappedTet.sourceTetIndex < sourceTetMesh->GetTetCount() &&
                sourceTetMesh->IsTetActive(mappedTet.sourceTetIndex);

            int destinationNodes[4];
            if (!TryGetTetNodes(
                    *destinationTetMesh,
                    destinationTetIndex,
                    destinationNodes))
            {
                destinationTetMesh->SetTetActive(destinationTetIndex, false);
                continue;
            }

            for (int destinationNode : destinationNodes)
            {
                if (!shouldBeActive)
                {
                    break;
                }

                const std::size_t vertexIndex =
                    static_cast<std::size_t>(destinationNode);
                if (vertexIndex >= embeddedDestinationVertices.size())
                {
                    shouldBeActive = false;
                    break;
                }

                const TetMeshMappedVertex& mappedVertex =
                    embeddedDestinationVertices[vertexIndex];
                shouldBeActive =
                    mappedVertex.valid &&
                    mappedVertex.sourceTetIndex >= 0 &&
                    mappedVertex.sourceTetIndex < sourceTetMesh->GetTetCount() &&
                    sourceTetMesh->IsTetActive(mappedVertex.sourceTetIndex);
            }

            destinationTetMesh->SetTetActive(destinationTetIndex, shouldBeActive);
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

            int sourceNodes[4];
            const bool sourceNodesAreValid =
                TryGetTetNodes(*sourceTetMesh, mappedVertex.sourceTetIndex, sourceNodes);

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
        activeStateDirty = true;
    }

    void TetMeshMapperComponent::MarkActiveStateDirty()
    {
        activeStateDirty = true;
    }

    bool TetMeshMapperComponent::IsMappingDirty() const
    {
        return mappingDirty;
    }

    void TetMeshMapperComponent::OnPhysicsEvent(const PhysicsEvent& event)
    {
        if (event.type != PhysicsEventType::TetMeshTopologyChanged)
        {
            return;
        }

        World* world = event.world;
        if (world == nullptr ||
            event.sender != world->GetComponent(sourceTetMeshHandle))
        {
            return;
        }

        MarkActiveStateDirty();
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
            activeStateDirty = true;
        }

        if (activeStateDirty)
        {
            if (!RefreshDestinationActiveStates(world))
            {
                return;
            }

            activeStateDirty = false;
        }

        UpdateDestinationNodes(world);
    }
}
