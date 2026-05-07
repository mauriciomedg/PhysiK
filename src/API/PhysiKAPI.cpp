#include "PhysiK/API/PhysiKAPI.h"

#include "PhysiK/Core/World/World.h"

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

    PHYSIK_API int PHYSIK_AddTet(
        PhysiK::WorldHandle world,
        int node0,
        int node1,
        int node2,
        int node3)
    {
        if (PhysiK::World* worldPtr = AsWorld(world))
        {
            return worldPtr->AddTet(node0, node1, node2, node3);
        }

        return -1;
    }

    PHYSIK_API PhysiK::ComponentHandle PHYSIK_CreateTetMeshComponent(
        PhysiK::WorldHandle world,
        const int* nodeIndices,
        int nodeCount,
        const int* tetIndices,
        int tetCount)
    {
        if (PhysiK::World* worldPtr = AsWorld(world))
        {
            return &worldPtr->CreateTetMeshComponent(nodeIndices, nodeCount, tetIndices, tetCount);
        }

        return nullptr;
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
}
