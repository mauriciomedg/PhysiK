#include "PhysiK/Components/VisualMeshComponent.h"

#include <cmath>
#include <cstddef>
#include <utility>

#include "PhysiK/Components/TetMeshComponent.h"
#include "PhysiK/Core/World/World.h"
#include "PhysiK/PhysicsData/Node.h"
#include "PhysiK/PhysicsData/Tet.h"

namespace PhysiK
{
    bool ComputeTetBarycentric(
        const Vec3& p,
        const Vec3& a,
        const Vec3& b,
        const Vec3& c,
        const Vec3& d,
        Vec4& outBarycentric)
    {
        const Vec3 ab = b - a;
        const Vec3 ac = c - a;
        const Vec3 ad = d - a;
        const Vec3 ap = p - a;
        const float determinant = Dot(ab, Cross(ac, ad));

        if (std::abs(determinant) <= 0.00000001f)
        {
            outBarycentric = Vec4{};
            return false;
        }

        const float inverseDeterminant = 1.0f / determinant;
        const float bWeight = Dot(ap, Cross(ac, ad)) * inverseDeterminant;
        const float cWeight = Dot(ab, Cross(ap, ad)) * inverseDeterminant;
        const float dWeight = Dot(ab, Cross(ac, ap)) * inverseDeterminant;
        const float aWeight = 1.0f - bWeight - cWeight - dWeight;

        outBarycentric = Vec4{aWeight, bWeight, cWeight, dWeight};
        return true;
    }

    VisualMeshComponent::VisualMeshComponent()
    {
        listenedEvents.push_back(PhysicsEventType::TetMeshTopologyChanged);
    }

    VisualMeshComponent::VisualMeshComponent(
        ComponentHandle hostTetMeshHandle,
        std::string debugEntityName)
        : VisualMeshComponent()
    {
        this->hostTetMeshHandle = hostTetMeshHandle;
        this->debugEntityName = std::move(debugEntityName);
    }

    void VisualMeshComponent::SetVisualMesh(
        const Vec3* vertices,
        int vertexCount,
        const int* triangleIndices,
        int triangleIndexCount)
    {
        restVisualVertices.clear();
        deformedVisualVertices.clear();
        this->triangleIndices.clear();
        filteredTriangleIndices.clear();
        embeddedVertices.clear();
        triangleValid.clear();

        if (vertices != nullptr && vertexCount > 0)
        {
            restVisualVertices.assign(vertices, vertices + vertexCount);
            deformedVisualVertices = restVisualVertices;
            embeddedVertices.resize(static_cast<std::size_t>(vertexCount));
        }

        if (triangleIndices != nullptr && triangleIndexCount > 0)
        {
            this->triangleIndices.assign(
                triangleIndices,
                triangleIndices + triangleIndexCount);
            filteredTriangleIndices = this->triangleIndices;
            triangleValid.assign(
                static_cast<std::size_t>(triangleIndexCount / 3),
                true);
        }
    }

    void VisualMeshComponent::BuildEmbedding(const World& world)
    {
        embeddedVertices.clear();
        embeddedVertices.resize(restVisualVertices.size());

        const Component* hostComponent = world.GetComponent(hostTetMeshHandle);
        const TetMeshComponent* hostTetMesh =
            dynamic_cast<const TetMeshComponent*>(hostComponent);
        if (hostTetMesh == nullptr)
        {
            return;
        }

        constexpr float insideEpsilon = -0.0001f;

        for (std::size_t vertexIndex = 0; vertexIndex < restVisualVertices.size(); ++vertexIndex)
        {
            EmbeddedVertex& embeddedVertex = embeddedVertices[vertexIndex];
            embeddedVertex = EmbeddedVertex{};

            for (int tetIndex = 0; tetIndex < hostTetMesh->GetTetCount(); ++tetIndex)
            {
                if (!hostTetMesh->IsTetActive(tetIndex))
                {
                    continue;
                }

                const int nodeIndices[4] = {
                    hostTetMesh->GetTetNodeIndex(tetIndex, 0),
                    hostTetMesh->GetTetNodeIndex(tetIndex, 1),
                    hostTetMesh->GetTetNodeIndex(tetIndex, 2),
                    hostTetMesh->GetTetNodeIndex(tetIndex, 3)};
                bool nodesAreValid = true;
                for (int nodeIndex : nodeIndices)
                {
                    if (nodeIndex < 0 || nodeIndex >= hostTetMesh->GetNodeCount())
                    {
                        nodesAreValid = false;
                    }
                }

                if (!nodesAreValid)
                {
                    continue;
                }

                Vec4 barycentric;
                if (!ComputeTetBarycentric(
                        restVisualVertices[vertexIndex],
                        hostTetMesh->GetLocalRestPosition(nodeIndices[0]),
                        hostTetMesh->GetLocalRestPosition(nodeIndices[1]),
                        hostTetMesh->GetLocalRestPosition(nodeIndices[2]),
                        hostTetMesh->GetLocalRestPosition(nodeIndices[3]),
                        barycentric))
                {
                    continue;
                }

                if (barycentric.x >= insideEpsilon &&
                    barycentric.y >= insideEpsilon &&
                    barycentric.z >= insideEpsilon &&
                    barycentric.w >= insideEpsilon)
                {
                    embeddedVertex.tetIndex = tetIndex;
                    embeddedVertex.barycentric = barycentric;
                    embeddedVertex.valid = true;
                    break;
                }
            }
        }

        UpdateTriangleValidity(world);
    }

    void VisualMeshComponent::UpdateTriangleValidity(const World& world)
    {
        filteredTriangleIndices.clear();

        const std::size_t triangleCount = triangleIndices.size() / 3u;
        triangleValid.assign(triangleCount, false);
        if (triangleCount == 0u)
        {
            return;
        }

        const Component* hostComponent = world.GetComponent(hostTetMeshHandle);
        const TetMeshComponent* hostTetMesh =
            dynamic_cast<const TetMeshComponent*>(hostComponent);
        if (hostTetMesh == nullptr)
        {
            return;
        }

        for (std::size_t triangleIndex = 0; triangleIndex < triangleCount; ++triangleIndex)
        {
            const std::size_t firstIndex = triangleIndex * 3u;
            const int vertex0 = triangleIndices[firstIndex];
            const int vertex1 = triangleIndices[firstIndex + 1u];
            const int vertex2 = triangleIndices[firstIndex + 2u];
            const int vertices[3] = {vertex0, vertex1, vertex2};

            bool valid = true;
            for (int vertexIndex : vertices)
            {
                if (vertexIndex < 0 ||
                    vertexIndex >= static_cast<int>(embeddedVertices.size()))
                {
                    valid = false;
                    break;
                }

                const EmbeddedVertex& embeddedVertex =
                    embeddedVertices[static_cast<std::size_t>(vertexIndex)];
                if (!embeddedVertex.valid ||
                    embeddedVertex.tetIndex < 0 ||
                    embeddedVertex.tetIndex >= hostTetMesh->GetTetCount() ||
                    !hostTetMesh->IsTetActive(embeddedVertex.tetIndex))
                {
                    valid = false;
                    break;
                }
            }

            triangleValid[triangleIndex] = valid;
            if (valid)
            {
                filteredTriangleIndices.push_back(vertex0);
                filteredTriangleIndices.push_back(vertex1);
                filteredTriangleIndices.push_back(vertex2);
            }
        }
    }

    void VisualMeshComponent::UpdateDeformedVertices(const World& world)
    {
        if (deformedVisualVertices.size() != restVisualVertices.size())
        {
            deformedVisualVertices = restVisualVertices;
        }

        if (embeddedVertices.size() != restVisualVertices.size())
        {
            return;
        }

        const Component* hostComponent = world.GetComponent(hostTetMeshHandle);
        const TetMeshComponent* hostTetMesh =
            dynamic_cast<const TetMeshComponent*>(hostComponent);
        if (hostTetMesh == nullptr)
        {
            return;
        }

        for (std::size_t vertexIndex = 0; vertexIndex < embeddedVertices.size(); ++vertexIndex)
        {
            const EmbeddedVertex& embeddedVertex = embeddedVertices[vertexIndex];
            if (!embeddedVertex.valid ||
                embeddedVertex.tetIndex < 0 ||
                embeddedVertex.tetIndex >= hostTetMesh->GetTetCount())
            {
                continue;
            }

            if (!hostTetMesh->IsTetActive(embeddedVertex.tetIndex))
            {
                continue;
            }

            const int nodeIndices[4] = {
                hostTetMesh->GetTetNodeIndex(embeddedVertex.tetIndex, 0),
                hostTetMesh->GetTetNodeIndex(embeddedVertex.tetIndex, 1),
                hostTetMesh->GetTetNodeIndex(embeddedVertex.tetIndex, 2),
                hostTetMesh->GetTetNodeIndex(embeddedVertex.tetIndex, 3)};
            bool nodesAreValid = true;
            for (int nodeIndex : nodeIndices)
            {
                if (nodeIndex < 0 || nodeIndex >= hostTetMesh->GetNodeCount())
                {
                    nodesAreValid = false;
                }
            }

            if (!nodesAreValid)
            {
                continue;
            }

            const Vec4& weights = embeddedVertex.barycentric;
            const int worldNode0 = hostTetMesh->GetWorldNodeIndex(nodeIndices[0]);
            const int worldNode1 = hostTetMesh->GetWorldNodeIndex(nodeIndices[1]);
            const int worldNode2 = hostTetMesh->GetWorldNodeIndex(nodeIndices[2]);
            const int worldNode3 = hostTetMesh->GetWorldNodeIndex(nodeIndices[3]);
            const Vec3& current0 = worldNode0 >= 0 ?
                world.GetNode(worldNode0).position :
                hostTetMesh->GetLocalCurrentPosition(nodeIndices[0]);
            const Vec3& current1 = worldNode1 >= 0 ?
                world.GetNode(worldNode1).position :
                hostTetMesh->GetLocalCurrentPosition(nodeIndices[1]);
            const Vec3& current2 = worldNode2 >= 0 ?
                world.GetNode(worldNode2).position :
                hostTetMesh->GetLocalCurrentPosition(nodeIndices[2]);
            const Vec3& current3 = worldNode3 >= 0 ?
                world.GetNode(worldNode3).position :
                hostTetMesh->GetLocalCurrentPosition(nodeIndices[3]);
            deformedVisualVertices[vertexIndex] =
                current0 * weights.x +
                current1 * weights.y +
                current2 * weights.z +
                current3 * weights.w;
        }
    }

    const std::vector<Vec3>& VisualMeshComponent::GetDeformedVertices() const
    {
        return deformedVisualVertices;
    }

    const std::vector<int>& VisualMeshComponent::GetTriangleIndices() const
    {
        return filteredTriangleIndices;
    }

    void VisualMeshComponent::OnPhysicsEvent(const PhysicsEvent& event)
    {
        if (event.type != PhysicsEventType::TetMeshTopologyChanged ||
            event.world == nullptr)
        {
            return;
        }

        if (event.world->GetComponent(hostTetMeshHandle) != event.sender)
        {
            return;
        }

        topologyDirty = true;
    }

    void VisualMeshComponent::PostUpdate(World& world, float dt)
    {
        (void)dt;

        if (topologyDirty)
        {
            UpdateTriangleValidity(world);
            topologyDirty = false;
        }

        UpdateDeformedVertices(world);
    }
}
