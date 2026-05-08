#include "PhysiK/API/PhysiKAPI.h"

#include "PhysiK/Components/CollisionComponent.h"
#include "PhysiK/Components/CollisionSphereComponent.h"
#include "PhysiK/Components/TetMeshComponent.h"
#include "PhysiK/Core/World/World.h"

#include <utility>

namespace
{
    PhysiK::World* AsWorld(PhysiK::WorldHandle handle)
    {
        return static_cast<PhysiK::World*>(handle);
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
        float z,
        float inverseMass)
    {
        if (PhysiK::World* worldPtr = AsWorld(world))
        {
            return worldPtr->AddNode(PhysiK::Vec3{x, y, z}, inverseMass);
        }

        return -1;
    }

    PHYSIK_API PhysiK::ComponentHandle PHYSIK_CreateTetMeshComponent(
        PhysiK::WorldHandle world,
        const int* nodeIndices,
        int nodeCount,
        const int* tetNodeIndices,
        int tetCount)
    {
        if (PhysiK::World* worldPtr = AsWorld(world))
        {
            PhysiK::Material material;
            auto component = PhysiK::TetMeshComponent::CreateFromGlobalNodes(
                *worldPtr,
                nodeIndices,
                nodeCount,
                tetNodeIndices,
                tetCount,
                material);
            return worldPtr->AddComponent(std::move(component));
        }

        return PhysiK::ComponentHandle{};
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

        return static_cast<int>(tetMesh->tets.size());
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
}
