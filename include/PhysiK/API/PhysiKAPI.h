#pragma once

#include "PhysiK/API/Handles.h"

#ifdef _WIN32
#if defined(PHYSIK_BUILDING_DLL)
#define PHYSIK_API __declspec(dllexport)
#else
#define PHYSIK_API __declspec(dllimport)
#endif
#else
#define PHYSIK_API
#endif

struct PhysikMaterialDesc
{
    float density;
    float youngModulus;
    float poissonRatio;
    float damping;
};

enum PhysikOverlapGeometryType
{
    PHYSIK_OverlapGeometry_Unknown = 0,
    PHYSIK_OverlapGeometry_Tetrahedron = 1,
    PHYSIK_OverlapGeometry_Triangle = 2,
    PHYSIK_OverlapGeometry_Sphere = 3,
    PHYSIK_OverlapGeometry_Node = 4
};

struct PhysikCollisionSphereOverlap
{
    int geometryType;

    PhysiK::ComponentHandle component;
    int primitiveIndex;

    int node0;
    int node1;
    int node2;
    int node3;

    int overlappedNodeMask;
    int overlappedNodeCount;

    float sphereCenterX;
    float sphereCenterY;
    float sphereCenterZ;
    float sphereRadius;
    float minDistance;
};

extern "C"
{
    PHYSIK_API PhysiK::WorldHandle PHYSIK_CreateWorld();
    PHYSIK_API void PHYSIK_DestroyWorld(PhysiK::WorldHandle world);

    PHYSIK_API void PHYSIK_Step(PhysiK::WorldHandle world, float dt);
    PHYSIK_API void PHYSIK_SetSubstepCount(PhysiK::WorldHandle world, int substepCount);
    PHYSIK_API void PHYSIK_SetSolverMode(PhysiK::WorldHandle world, int mode);
    PHYSIK_API void PHYSIK_SetGravity(
        PhysiK::WorldHandle world,
        float x,
        float y,
        float z);
    PHYSIK_API void PHYSIK_SetExternalLogicCallback(
        PhysiK::WorldHandle world,
        PhysiK::ExternalLogicCallback callback,
        void* userData);
    PHYSIK_API void PHYSIK_ClearExternalLogicCallback(PhysiK::WorldHandle world);

    PHYSIK_API int PHYSIK_AddNode(
        PhysiK::WorldHandle world,
        float x,
        float y,
        float z);

    PHYSIK_API void PHYSIK_SetNodeFixed(
        PhysiK::WorldHandle world,
        int nodeIndex,
        int fixed);

    PHYSIK_API int PHYSIK_IsNodeFixed(
        PhysiK::WorldHandle world,
        int nodeIndex);

    PHYSIK_API PhysiK::ComponentHandle PHYSIK_CreateTetMeshComponent(
        PhysiK::WorldHandle world,
        const int* nodeIndices,
        int nodeCount,
        const int* tetNodeIndices,
        int tetCount,
        const PhysikMaterialDesc* material,
        int femModel);

    PHYSIK_API void PHYSIK_SetTetMeshMaterial(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle component,
        const PhysikMaterialDesc* material);

    PHYSIK_API PhysiK::ComponentHandle PHYSIK_CreateCollisionSphereComponent(
        PhysiK::WorldHandle world,
        float x,
        float y,
        float z,
        float radius);

    PHYSIK_API void PHYSIK_DestroyComponent(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle component);

    PHYSIK_API int PHYSIK_IsComponentHandleValid(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle component);

    PHYSIK_API int PHYSIK_GetTetMeshTetCount(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle component);

    PHYSIK_API int PHYSIK_IsTetActive(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle component,
        int tetIndex);

    PHYSIK_API void PHYSIK_SetTetActive(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle component,
        int tetIndex,
        int active);

    PHYSIK_API void PHYSIK_DeactivateTet(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle component,
        int tetIndex);

    PHYSIK_API int PHYSIK_GetActiveTetCount(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle component);

    PHYSIK_API int PHYSIK_GetCollisionSphereOverlapCount(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle sphereComponent);

    PHYSIK_API int PHYSIK_GetCollisionSphereOverlaps(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle sphereComponent,
        PhysikCollisionSphereOverlap* outOverlaps,
        int maxOverlaps);

    PHYSIK_API void PHYSIK_SetCollisionComponentKinematicTarget(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle component,
        float x,
        float y,
        float z);

    PHYSIK_API int PHYSIK_GetPointConnectionCount(PhysiK::WorldHandle world);

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
        float damping);

    PHYSIK_API void PHYSIK_GetNodePosition(
        PhysiK::WorldHandle world,
        int nodeIndex,
        float* outX,
        float* outY,
        float* outZ);

    PHYSIK_API void PHYSIK_GetNodeVelocity(
        PhysiK::WorldHandle world,
        int nodeIndex,
        float* outX,
        float* outY,
        float* outZ);

    PHYSIK_API void PHYSIK_SetNodePosition(
        PhysiK::WorldHandle world,
        int nodeIndex,
        float x,
        float y,
        float z);

    PHYSIK_API void PHYSIK_SetNodeVelocity(
        PhysiK::WorldHandle world,
        int nodeIndex,
        float x,
        float y,
        float z);
}
