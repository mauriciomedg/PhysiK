#include "PhysiK/API/PhysiKAPI.h"

#include "PhysiK/Components/CollisionComponent.h"
#include "PhysiK/Components/CollisionSphereComponent.h"
#include "PhysiK/Components/SurfaceExtractionComponent.h"
#include "PhysiK/Components/TetMeshComponent.h"
#include "PhysiK/Components/TetMeshPhysicsComponent.h"
#include "PhysiK/Components/TetMeshMapperComponent.h"
#include "PhysiK/Components/VisualMeshComponent.h"
#include "PhysiK/Core/World/World.h"

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

namespace
{
    PhysiK::World* AsWorld(PhysiK::WorldHandle handle)
    {
        return static_cast<PhysiK::World*>(handle);
    }

    PhysiK::Material ToMaterial(const PhysikMaterialDesc& desc)
    {
        PhysiK::Material material;
        material.density = desc.density;
        material.youngModulus = desc.youngModulus;
        material.poissonRatio = desc.poissonRatio;
        material.damping = desc.damping;
        return material;
    }

    PhysiK::FemModel ToFemModel(int value)
    {
        switch (value)
        {
        case 1:
            return PhysiK::FemModel::Corotational;
        case 2:
            return PhysiK::FemModel::NeoHookean;
        case 0:
        default:
            return PhysiK::FemModel::Linear;
        }
    }

    PhysikCollisionSphereOverlap ToApiOverlap(const PhysiK::CollisionSphereOverlap& overlap)
    {
        PhysikCollisionSphereOverlap apiOverlap;
        apiOverlap.geometryType = static_cast<int>(overlap.geometryType);
        apiOverlap.component = overlap.component;
        apiOverlap.primitiveIndex = overlap.primitiveIndex;
        apiOverlap.node0 = overlap.node0;
        apiOverlap.node1 = overlap.node1;
        apiOverlap.node2 = overlap.node2;
        apiOverlap.node3 = overlap.node3;
        apiOverlap.overlappedNodeMask = overlap.overlappedNodeMask;
        apiOverlap.overlappedNodeCount = overlap.overlappedNodeCount;
        apiOverlap.sphereCenterX = overlap.sphereCenter.x;
        apiOverlap.sphereCenterY = overlap.sphereCenter.y;
        apiOverlap.sphereCenterZ = overlap.sphereCenter.z;
        apiOverlap.sphereRadius = overlap.sphereRadius;
        apiOverlap.minDistance = overlap.minDistance;
        return apiOverlap;
    }

    std::vector<PhysiK::CollisionSphereOverlap> QueryCollisionSphereOverlaps(
        PhysiK::World* world,
        PhysiK::ComponentHandle sphereComponent)
    {
        std::vector<PhysiK::CollisionSphereOverlap> overlaps;
        if (world == nullptr)
        {
            return overlaps;
        }

        const auto* sphere = dynamic_cast<const PhysiK::CollisionSphereComponent*>(
            world->GetComponent(sphereComponent));
        if (sphere == nullptr)
        {
            return overlaps;
        }

        sphere->QueryOverlaps(*world, overlaps);
        return overlaps;
    }

    PhysiK::CollisionSphereComponent* AsCollisionSphere(
        PhysiK::World* world,
        PhysiK::ComponentHandle sphereComponent)
    {
        if (world == nullptr)
        {
            return nullptr;
        }

        return dynamic_cast<PhysiK::CollisionSphereComponent*>(
            world->GetComponent(sphereComponent));
    }

    PhysiK::TetMeshComponent* AsTetMesh(
        PhysiK::World* world,
        PhysiK::ComponentHandle tetMeshComponent)
    {
        if (world == nullptr)
        {
            return nullptr;
        }

        return dynamic_cast<PhysiK::TetMeshComponent*>(
            world->GetComponent(tetMeshComponent));
    }

    PhysiK::TetMeshPhysicsComponent* AsTetMeshPhysics(
        PhysiK::World* world,
        PhysiK::ComponentHandle tetMeshComponent)
    {
        if (world == nullptr)
        {
            return nullptr;
        }

        return dynamic_cast<PhysiK::TetMeshPhysicsComponent*>(
            world->GetComponent(tetMeshComponent));
    }

    PhysiK::VisualMeshComponent* AsVisualMesh(
        PhysiK::World* world,
        PhysiK::ComponentHandle visualMeshComponent)
    {
        if (world == nullptr)
        {
            return nullptr;
        }

        return dynamic_cast<PhysiK::VisualMeshComponent*>(
            world->GetComponent(visualMeshComponent));
    }

    const PhysiK::VisualMeshComponent* AsVisualMesh(
        const PhysiK::World* world,
        PhysiK::ComponentHandle visualMeshComponent)
    {
        if (world == nullptr)
        {
            return nullptr;
        }

        return dynamic_cast<const PhysiK::VisualMeshComponent*>(
            world->GetComponent(visualMeshComponent));
    }

    PhysiK::TetMeshMapperComponent* AsTetMeshMapper(
        PhysiK::World* world,
        PhysiK::ComponentHandle mapperComponent)
    {
        if (world == nullptr)
        {
            return nullptr;
        }

        return dynamic_cast<PhysiK::TetMeshMapperComponent*>(
            world->GetComponent(mapperComponent));
    }
}

extern "C"
{
    PHYSIK_API PhysiK::WorldHandle PHYSIK_CreateWorld()
    {
        return new PhysiK::World();
    }

    PHYSIK_API void PHYSIK_DestroyWorld(PhysiK::WorldHandle world)
    {
        delete AsWorld(world);
    }

    PHYSIK_API void PHYSIK_Step(PhysiK::WorldHandle world, float dt)
    {
        if (PhysiK::World* worldPtr = AsWorld(world))
        {
            worldPtr->Step(dt);
        }
    }

    PHYSIK_API void PHYSIK_SetSubstepCount(PhysiK::WorldHandle world, int substepCount)
    {
        if (PhysiK::World* worldPtr = AsWorld(world))
        {
            worldPtr->SetSubstepCount(substepCount);
        }
    }

    PHYSIK_API void PHYSIK_SetSolverMode(PhysiK::WorldHandle world, int mode)
    {
        if (PhysiK::World* worldPtr = AsWorld(world))
        {
            worldPtr->SetSolverMode(
                mode == 1 ? PhysiK::SolverMode::ImplicitEuler : PhysiK::SolverMode::Explicit);
        }
    }

    PHYSIK_API void PHYSIK_EnablePerformanceLogging(
        PhysiK::WorldHandle world,
        int enabled)
    {
        if (PhysiK::World* worldPtr = AsWorld(world))
        {
            worldPtr->EnablePerformanceLogging(enabled != 0);
        }
    }

    PHYSIK_API void PHYSIK_SetPerformanceLogPath(
        PhysiK::WorldHandle world,
        const char* path)
    {
        if (PhysiK::World* worldPtr = AsWorld(world))
        {
            worldPtr->SetPerformanceLogPath(path);
        }
    }

    PHYSIK_API void PHYSIK_SetGravity(
        PhysiK::WorldHandle world,
        float x,
        float y,
        float z)
    {
        if (PhysiK::World* worldPtr = AsWorld(world))
        {
            worldPtr->SetGravity(PhysiK::Vec3{x, y, z});
        }
    }

    PHYSIK_API void PHYSIK_SetExternalLogicCallback(
        PhysiK::WorldHandle world,
        PhysiK::ExternalLogicCallback callback,
        void* userData)
    {
        if (PhysiK::World* worldPtr = AsWorld(world))
        {
            worldPtr->SetExternalLogicCallback(callback, userData);
        }
    }

    PHYSIK_API void PHYSIK_ClearExternalLogicCallback(PhysiK::WorldHandle world)
    {
        if (PhysiK::World* worldPtr = AsWorld(world))
        {
            worldPtr->ClearExternalLogicCallback();
        }
    }

    PHYSIK_API int PHYSIK_AddNode(
        PhysiK::WorldHandle world,
        float x,
        float y,
        float z)
    {
        if (PhysiK::World* worldPtr = AsWorld(world))
        {
            return worldPtr->AddNode(PhysiK::Vec3{x, y, z});
        }

        return -1;
    }

    PHYSIK_API void PHYSIK_SetNodeFixed(
        PhysiK::WorldHandle world,
        int nodeIndex,
        int fixed)
    {
        if (PhysiK::World* worldPtr = AsWorld(world))
        {
            worldPtr->SetNodeFixed(nodeIndex, fixed != 0);
        }
    }

    PHYSIK_API int PHYSIK_IsNodeFixed(
        PhysiK::WorldHandle world,
        int nodeIndex)
    {
        if (PhysiK::World* worldPtr = AsWorld(world))
        {
            return worldPtr->IsNodeFixed(nodeIndex) ? 1 : 0;
        }

        return 0;
    }

    PHYSIK_API PhysiK::ComponentHandle PHYSIK_CreateTetMeshComponent(
        PhysiK::WorldHandle world,
        const PhysiK::Vec3* positions,
        int nodeCount,
        const int* tetLocalNodeIndices,
        int tetCount)
    {
        PhysiK::World* worldPtr = AsWorld(world);
        if (worldPtr == nullptr ||
            positions == nullptr ||
            nodeCount <= 0 ||
            tetLocalNodeIndices == nullptr ||
            tetCount <= 0)
        {
            return PhysiK::ComponentHandle{};
        }

        auto component = std::make_unique<PhysiK::TetMeshComponent>();

        component->SetGeometry(
            positions,
            nodeCount,
            tetLocalNodeIndices,
            tetCount);

        return worldPtr->AddComponent(std::move(component));
    }

    PHYSIK_API PhysiK::ComponentHandle PHYSIK_CreateTetMeshPhysicsComponent(
        PhysiK::WorldHandle world,
        const int* nodeIndices,
        int nodeCount,
        const int* tetNodeIndices,
        int tetCount,
        const PhysikMaterialDesc* material,
        int femModel)
    {
        PhysiK::World* worldPtr = AsWorld(world);
        if (worldPtr == nullptr || material == nullptr)
        {
            return PhysiK::ComponentHandle{};
        }

        PhysiK::TetMeshPhysicsComponentDesc desc;
        desc.material = ToMaterial(*material);
        desc.femModel = ToFemModel(femModel);

        auto component = PhysiK::TetMeshPhysicsComponent::CreateFromGlobalNodes(
            *worldPtr,
            nodeIndices,
            nodeCount,
            tetNodeIndices,
            tetCount,
            desc);
        return worldPtr->AddComponent(std::move(component));
    }

    PHYSIK_API PhysiK::ComponentHandle PHYSIK_CreateCollisionSphereComponent(
        PhysiK::WorldHandle world,
        float x,
        float y,
        float z,
        float radius)
    {
        if (PhysiK::World* worldPtr = AsWorld(world))
        {
            auto component = PhysiK::CollisionSphereComponent::Create(
                PhysiK::Vec3{x, y, z},
                radius);
            return worldPtr->AddComponent(std::move(component));
        }

        return PhysiK::ComponentHandle{};
    }

    PHYSIK_API PhysiK::ComponentHandle PHYSIK_CreateVisualMeshComponent(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle hostTetMesh)
    {
        PhysiK::World* worldPtr = AsWorld(world);
        if (worldPtr == nullptr || AsTetMesh(worldPtr, hostTetMesh) == nullptr)
        {
            return PhysiK::ComponentHandle{};
        }

        auto component = std::make_unique<PhysiK::VisualMeshComponent>(
            hostTetMesh,
            std::string{});
        return worldPtr->AddComponent(std::move(component));
    }

    PHYSIK_API PhysiK::ComponentHandle PHYSIK_CreateSurfaceExtractionComponent(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle hostTetMesh)
    {
        PhysiK::World* worldPtr = AsWorld(world);
        if (worldPtr == nullptr || AsTetMesh(worldPtr, hostTetMesh) == nullptr)
        {
            return PhysiK::ComponentHandle{};
        }

        auto component =
            std::make_unique<PhysiK::SurfaceExtractionComponent>(hostTetMesh);
        return worldPtr->AddComponent(std::move(component));
    }

    PHYSIK_API PhysiK::ComponentHandle PHYSIK_CreateTetMeshMapperComponent(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle sourceTetMesh,
        PhysiK::ComponentHandle destinationTetMesh)
    {
        PhysiK::World* worldPtr = AsWorld(world);
        if (worldPtr == nullptr ||
            AsTetMesh(worldPtr, sourceTetMesh) == nullptr ||
            AsTetMesh(worldPtr, destinationTetMesh) == nullptr)
        {
            return PhysiK::ComponentHandle{};
        }

        auto component = std::make_unique<PhysiK::TetMeshMapperComponent>(
            sourceTetMesh,
            destinationTetMesh);
        return worldPtr->AddComponent(std::move(component));
    }

    PHYSIK_API int PHYSIK_BuildTetMeshMapping(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle mapper)
    {
        PhysiK::World* worldPtr = AsWorld(world);
        PhysiK::TetMeshMapperComponent* mapperComponent =
            AsTetMeshMapper(worldPtr, mapper);
        if (worldPtr == nullptr || mapperComponent == nullptr)
        {
            return 0;
        }

        mapperComponent->MarkMappingDirty();
        return 1;
    }

    PHYSIK_API int PHYSIK_SetTetMeshLocalCurrentPosition(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle tetMesh,
        int localNodeIndex,
        float x,
        float y,
        float z)
    {
        PhysiK::World* worldPtr = AsWorld(world);
        if (worldPtr == nullptr)
        {
            return 0;
        }

        auto* tetMeshComponent = dynamic_cast<PhysiK::TetMeshComponent*>(
            worldPtr->GetComponent(tetMesh));
        if (tetMeshComponent == nullptr ||
            localNodeIndex < 0 ||
            localNodeIndex >= tetMeshComponent->GetNodeCount())
        {
            return 0;
        }

        tetMeshComponent->SetLocalCurrentPosition(
            localNodeIndex,
            PhysiK::Vec3{x, y, z});
        return 1;
    }

    PHYSIK_API int PHYSIK_GetTetMeshLocalCurrentPosition(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle tetMesh,
        int localNodeIndex,
        float* outX,
        float* outY,
        float* outZ)
    {
        PhysiK::World* worldPtr = AsWorld(world);
        if (worldPtr == nullptr ||
            outX == nullptr ||
            outY == nullptr ||
            outZ == nullptr)
        {
            return 0;
        }

        const auto* tetMeshComponent = dynamic_cast<const PhysiK::TetMeshComponent*>(
            worldPtr->GetComponent(tetMesh));
        if (tetMeshComponent == nullptr ||
            localNodeIndex < 0 ||
            localNodeIndex >= tetMeshComponent->GetNodeCount())
        {
            return 0;
        }

        const PhysiK::Vec3& position =
            tetMeshComponent->GetLocalCurrentPosition(localNodeIndex);
        *outX = position.x;
        *outY = position.y;
        *outZ = position.z;
        return 1;
    }

    PHYSIK_API void PHYSIK_SetVisualMeshData(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle visualMesh,
        const PhysiK::Vec3* vertices,
        int vertexCount,
        const int* triangleIndices,
        int triangleIndexCount)
    {
        PhysiK::VisualMeshComponent* visual =
            AsVisualMesh(AsWorld(world), visualMesh);
        if (visual == nullptr)
        {
            return;
        }

        visual->SetVisualMesh(vertices, vertexCount, triangleIndices, triangleIndexCount);
    }

    PHYSIK_API int PHYSIK_BuildVisualMeshEmbedding(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle visualMesh)
    {
        PhysiK::World* worldPtr = AsWorld(world);
        PhysiK::VisualMeshComponent* visual = AsVisualMesh(worldPtr, visualMesh);
        if (worldPtr == nullptr || visual == nullptr)
        {
            return 0;
        }

        visual->BuildEmbedding(*worldPtr);
        for (const PhysiK::EmbeddedVertex& embeddedVertex : visual->embeddedVertices)
        {
            if (embeddedVertex.valid)
            {
                return 1;
            }
        }

        return 0;
    }

    PHYSIK_API int PHYSIK_GetVisualMeshVertexCount(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle visualMesh)
    {
        const PhysiK::VisualMeshComponent* visual =
            AsVisualMesh(AsWorld(world), visualMesh);
        if (visual == nullptr)
        {
            return 0;
        }

        return static_cast<int>(visual->GetDeformedVertices().size());
    }

    PHYSIK_API int PHYSIK_GetVisualMeshTriangleIndexCount(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle visualMesh)
    {
        const PhysiK::VisualMeshComponent* visual =
            AsVisualMesh(AsWorld(world), visualMesh);
        if (visual == nullptr)
        {
            return 0;
        }

        return static_cast<int>(visual->GetTriangleIndices().size());
    }

    PHYSIK_API int PHYSIK_CopyVisualMeshVertices(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle visualMesh,
        PhysiK::Vec3* outVertices,
        int maxVertexCount)
    {
        if (outVertices == nullptr || maxVertexCount <= 0)
        {
            return 0;
        }

        const PhysiK::VisualMeshComponent* visual =
            AsVisualMesh(AsWorld(world), visualMesh);
        if (visual == nullptr)
        {
            return 0;
        }

        const std::vector<PhysiK::Vec3>& vertices = visual->GetDeformedVertices();
        const int writeCount = std::min(maxVertexCount, static_cast<int>(vertices.size()));
        for (int i = 0; i < writeCount; ++i)
        {
            outVertices[i] = vertices[static_cast<std::size_t>(i)];
        }

        return writeCount;
    }

    PHYSIK_API int PHYSIK_CopyVisualMeshTriangleIndices(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle visualMesh,
        int* outIndices,
        int maxIndexCount)
    {
        if (outIndices == nullptr || maxIndexCount <= 0)
        {
            return 0;
        }

        const PhysiK::VisualMeshComponent* visual =
            AsVisualMesh(AsWorld(world), visualMesh);
        if (visual == nullptr)
        {
            return 0;
        }

        const std::vector<int>& indices = visual->GetTriangleIndices();
        const int writeCount = std::min(maxIndexCount, static_cast<int>(indices.size()));
        for (int i = 0; i < writeCount; ++i)
        {
            outIndices[i] = indices[static_cast<std::size_t>(i)];
        }

        return writeCount;
    }

    PHYSIK_API void PHYSIK_SetCollisionSphereConnectionSettings(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle sphereComponent,
        float stiffness,
        float damping)
    {
        PhysiK::CollisionSphereComponent* sphere =
            AsCollisionSphere(AsWorld(world), sphereComponent);
        if (sphere == nullptr)
        {
            return;
        }

        sphere->SetConnectionSettings(stiffness, damping);
    }

    PHYSIK_API void PHYSIK_GetCollisionSphereConnectionSettings(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle sphereComponent,
        float* outStiffness,
        float* outDamping)
    {
        const PhysiK::CollisionSphereComponent* sphere =
            AsCollisionSphere(AsWorld(world), sphereComponent);
        if (sphere == nullptr)
        {
            return;
        }

        if (outStiffness != nullptr)
        {
            *outStiffness = sphere->GetConnectionStiffness();
        }

        if (outDamping != nullptr)
        {
            *outDamping = sphere->GetConnectionDamping();
        }
    }

    PHYSIK_API void PHYSIK_SetTetMeshMaterial(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle component,
        const PhysikMaterialDesc* material)
    {
        PhysiK::World* worldPtr = AsWorld(world);
        if (worldPtr == nullptr || material == nullptr)
        {
            return;
        }

        auto* tetMesh = dynamic_cast<PhysiK::TetMeshPhysicsComponent*>(
            worldPtr->GetComponent(component));
        if (tetMesh == nullptr)
        {
            return;
        }

        tetMesh->SetMaterial(ToMaterial(*material));
    }

    PHYSIK_API void PHYSIK_DestroyComponent(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle component)
    {
        if (PhysiK::World* worldPtr = AsWorld(world))
        {
            worldPtr->DestroyComponent(component);
        }
    }

    PHYSIK_API int PHYSIK_IsComponentHandleValid(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle component)
    {
        if (PhysiK::World* worldPtr = AsWorld(world))
        {
            return worldPtr->GetComponent(component) != nullptr ? 1 : 0;
        }

        return 0;
    }

    PHYSIK_API int PHYSIK_GetTetMeshTetCount(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle component)
    {
        PhysiK::World* worldPtr = AsWorld(world);
        if (worldPtr == nullptr)
        {
            return 0;
        }

        const auto* tetMesh = dynamic_cast<const PhysiK::TetMeshComponent*>(
            worldPtr->GetComponent(component));
        if (tetMesh == nullptr)
        {
            return 0;
        }

        return tetMesh->GetTetCount();
    }

    PHYSIK_API int PHYSIK_IsTetActive(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle component,
        int tetIndex)
    {
        PhysiK::World* worldPtr = AsWorld(world);
        if (worldPtr == nullptr)
        {
            return 0;
        }

        const auto* tetMesh = dynamic_cast<const PhysiK::TetMeshComponent*>(
            worldPtr->GetComponent(component));
        if (tetMesh == nullptr)
        {
            return 0;
        }

        return tetMesh->IsTetActive(tetIndex) ? 1 : 0;
    }

    PHYSIK_API void PHYSIK_SetTetActive(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle component,
        int tetIndex,
        int active)
    {
        PhysiK::World* worldPtr = AsWorld(world);
        if (worldPtr == nullptr)
        {
            return;
        }

        auto* tetMesh = dynamic_cast<PhysiK::TetMeshComponent*>(worldPtr->GetComponent(component));
        if (tetMesh == nullptr)
        {
            return;
        }

        tetMesh->SetTetActive(tetIndex, active != 0);
    }

    PHYSIK_API void PHYSIK_DeactivateTet(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle component,
        int tetIndex)
    {
        PHYSIK_SetTetActive(world, component, tetIndex, 0);
    }

    PHYSIK_API int PHYSIK_GetActiveTetCount(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle component)
    {
        PhysiK::World* worldPtr = AsWorld(world);
        if (worldPtr == nullptr)
        {
            return 0;
        }

        const auto* tetMesh = dynamic_cast<const PhysiK::TetMeshComponent*>(
            worldPtr->GetComponent(component));
        if (tetMesh == nullptr)
        {
            return 0;
        }

        return tetMesh->GetActiveTetCount();
    }

    PHYSIK_API int PHYSIK_GetCollisionSphereOverlapCount(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle sphereComponent)
    {
        const std::vector<PhysiK::CollisionSphereOverlap> overlaps =
            QueryCollisionSphereOverlaps(AsWorld(world), sphereComponent);
        return static_cast<int>(overlaps.size());
    }

    PHYSIK_API int PHYSIK_GetCollisionSphereOverlaps(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle sphereComponent,
        PhysikCollisionSphereOverlap* outOverlaps,
        int maxOverlaps)
    {
        if (outOverlaps == nullptr || maxOverlaps <= 0)
        {
            return 0;
        }

        const std::vector<PhysiK::CollisionSphereOverlap> overlaps =
            QueryCollisionSphereOverlaps(AsWorld(world), sphereComponent);
        const int writeCount = std::min(maxOverlaps, static_cast<int>(overlaps.size()));
        for (int i = 0; i < writeCount; ++i)
        {
            outOverlaps[i] = ToApiOverlap(overlaps[static_cast<std::size_t>(i)]);
        }

        return writeCount;
    }

    PHYSIK_API void PHYSIK_SetCollisionComponentKinematicTarget(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle component,
        float x,
        float y,
        float z)
    {
        PhysiK::World* worldPtr = AsWorld(world);
        if (worldPtr == nullptr)
        {
            return;
        }

        auto* collision = dynamic_cast<PhysiK::CollisionComponent*>(worldPtr->GetComponent(component));
        if (collision == nullptr)
        {
            return;
        }

        PhysiK::Transform target = collision->transform;
        target.position = PhysiK::Vec3{x, y, z};
        collision->SetKinematicTarget(target);
    }

    PHYSIK_API int PHYSIK_GetPointConnectionCount(PhysiK::WorldHandle world)
    {
        if (PhysiK::World* worldPtr = AsWorld(world))
        {
            return worldPtr->GetTransientConnectionCount();
        }

        return 0;
    }

    PHYSIK_API void PHYSIK_AddPointConnection(
        PhysiK::WorldHandle world,
        int node0,
        int node1,
        int node2,
        int node3,
        float barycentricX,
        float barycentricY,
        float barycentricZ,
        float barycentricW,
        float targetX,
        float targetY,
        float targetZ,
        float stiffness,
        float damping)
    {
        if (PhysiK::World* worldPtr = AsWorld(world))
        {
            PhysiK::PointConnection connection;
            connection.node0 = node0;
            connection.node1 = node1;
            connection.node2 = node2;
            connection.node3 = node3;
            connection.barycentric = PhysiK::Vec4{
                barycentricX,
                barycentricY,
                barycentricZ,
                barycentricW};
            connection.targetPosition = PhysiK::Vec3{targetX, targetY, targetZ};
            connection.stiffness = stiffness;
            connection.damping = damping;

            worldPtr->AddPointConnection(connection);
        }
    }

    PHYSIK_API void PHYSIK_GetNodePosition(
        PhysiK::WorldHandle world,
        int nodeIndex,
        float* outX,
        float* outY,
        float* outZ)
    {
        if (PhysiK::World* worldPtr = AsWorld(world))
        {
            const PhysiK::Node& node = worldPtr->GetNode(nodeIndex);

            if (outX != nullptr)
            {
                *outX = node.position.x;
            }

            if (outY != nullptr)
            {
                *outY = node.position.y;
            }

            if (outZ != nullptr)
            {
                *outZ = node.position.z;
            }
        }
    }

    PHYSIK_API void PHYSIK_GetNodeVelocity(
        PhysiK::WorldHandle world,
        int nodeIndex,
        float* outX,
        float* outY,
        float* outZ)
    {
        if (PhysiK::World* worldPtr = AsWorld(world))
        {
            const PhysiK::Node& node = worldPtr->GetNode(nodeIndex);

            if (outX != nullptr)
            {
                *outX = node.velocity.x;
            }

            if (outY != nullptr)
            {
                *outY = node.velocity.y;
            }

            if (outZ != nullptr)
            {
                *outZ = node.velocity.z;
            }
        }
    }

    PHYSIK_API void PHYSIK_SetNodePosition(
        PhysiK::WorldHandle world,
        int nodeIndex,
        float x,
        float y,
        float z)
    {
        if (PhysiK::World* worldPtr = AsWorld(world))
        {
            worldPtr->SetNodePosition(nodeIndex, PhysiK::Vec3{x, y, z});
        }
    }

    PHYSIK_API void PHYSIK_SetNodeVelocity(
        PhysiK::WorldHandle world,
        int nodeIndex,
        float x,
        float y,
        float z)
    {
        if (PhysiK::World* worldPtr = AsWorld(world))
        {
            worldPtr->GetNode(nodeIndex).velocity = PhysiK::Vec3{x, y, z};
        }
    }
}
