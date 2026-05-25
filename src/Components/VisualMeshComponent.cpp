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

        const std::vector<Node>& nodes = world.GetNodes();
        constexpr float insideEpsilon = -0.0001f;

        for (std::size_t vertexIndex = 0; vertexIndex < restVisualVertices.size(); ++vertexIndex)
        {
            EmbeddedVertex& embeddedVertex = embeddedVertices[vertexIndex];
            embeddedVertex = EmbeddedVertex{};

            for (std::size_t tetIndex = 0; tetIndex < hostTetMesh->tets.size(); ++tetIndex)
            {
                const Tet& tet = hostTetMesh->tets[tetIndex];
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
                    embeddedVertex.tetIndex = static_cast<int>(tetIndex);
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
                    embeddedVertex.tetIndex >= static_cast<int>(hostTetMesh->tets.size()) ||
                    !hostTetMesh->tets[static_cast<std::size_t>(embeddedVertex.tetIndex)].active)
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

        const std::vector<Node>& nodes = world.GetNodes();

        for (std::size_t vertexIndex = 0; vertexIndex < embeddedVertices.size(); ++vertexIndex)
        {
            const EmbeddedVertex& embeddedVertex = embeddedVertices[vertexIndex];
            if (!embeddedVertex.valid ||
                embeddedVertex.tetIndex < 0 ||
                embeddedVertex.tetIndex >= static_cast<int>(hostTetMesh->tets.size()))
            {
                continue;
            }

            const Tet& tet =
                hostTetMesh->tets[static_cast<std::size_t>(embeddedVertex.tetIndex)];
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
