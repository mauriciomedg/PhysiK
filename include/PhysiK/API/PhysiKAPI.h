#pragma once

#include "PhysiK/API/Handles.h"
#include "PhysiK/Math/Vec3.h"

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
    PHYSIK_API void PHYSIK_SetConjugateGradientTolerance(
        PhysiK::WorldHandle world,
        float tolerance);
    PHYSIK_API float PHYSIK_GetConjugateGradientTolerance(
        PhysiK::WorldHandle world);
    PHYSIK_API void PHYSIK_SetConjugateGradientMaxIterations(
        PhysiK::WorldHandle world,
        int maxIterations);
    PHYSIK_API int PHYSIK_GetConjugateGradientMaxIterations(
        PhysiK::WorldHandle world);
    PHYSIK_API int PHYSIK_GetLastConjugateGradientIterations(
        PhysiK::WorldHandle world);
    PHYSIK_API float PHYSIK_GetLastConjugateGradientResidualNorm(
        PhysiK::WorldHandle world);
    PHYSIK_API int PHYSIK_DidLastConjugateGradientSolveConverge(
        PhysiK::WorldHandle world);
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

    PHYSIK_API PhysiK::GeneratedTetMeshHandle PHYSIK_GenerateTetMesh(
        const PhysiK::Vec3* positions,
        int nodeCount,
        const int* tetLocalNodeIndices,
        int tetCount);

    PHYSIK_API int PHYSIK_IsGeneratedTetMeshHandleValid(
        PhysiK::GeneratedTetMeshHandle generatedTetMeshHandle);

    PHYSIK_API void PHYSIK_DestroyGeneratedTetMesh(
        PhysiK::GeneratedTetMeshHandle generatedTetMeshHandle);

    PHYSIK_API int PHYSIK_GetGeneratedTetMeshVertexCount(
        PhysiK::GeneratedTetMeshHandle generatedTetMeshHandle);

    PHYSIK_API int PHYSIK_GetGeneratedTetMeshTetCount(
        PhysiK::GeneratedTetMeshHandle generatedTetMeshHandle);

    PHYSIK_API int PHYSIK_GetGeneratedTetMeshTetIndexCount(
        PhysiK::GeneratedTetMeshHandle generatedTetMeshHandle);

    PHYSIK_API int PHYSIK_GetGeneratedTetMeshVertex(
        PhysiK::GeneratedTetMeshHandle generatedTetMeshHandle,
        int vertexIndex,
        float* outX,
        float* outY,
        float* outZ);

    PHYSIK_API int PHYSIK_GetGeneratedTetMeshTetNodeIndex(
        PhysiK::GeneratedTetMeshHandle generatedTetMeshHandle,
        int tetIndexArrayIndex,
        int* outNodeIndex);

    PHYSIK_API PhysiK::ComponentHandle PHYSIK_CreateTetMeshComponent(
        PhysiK::WorldHandle world,
        PhysiK::GeneratedTetMeshHandle generatedTetMeshHandle);

    PHYSIK_API PhysiK::ComponentHandle PHYSIK_CreateTetMeshPhysicsComponent(
        PhysiK::WorldHandle world,
        PhysiK::GeneratedTetMeshHandle generatedTetMeshHandle,
        const PhysikMaterialDesc* material);

    PHYSIK_API int PHYSIK_GetTetMeshGlobalNodeBeginIndex(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle tetMeshPhysicsHandle);

    PHYSIK_API int PHYSIK_GetTetMeshGlobalNodeCount(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle tetMeshPhysicsHandle);

    PHYSIK_API int PHYSIK_GetTetMeshGlobalNodeIndex(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle tetMeshPhysicsHandle,
        int localNodeIndex);

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

    PHYSIK_API PhysiK::ComponentHandle PHYSIK_CreateSurfaceExtractionComponent(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle hostTetMesh);

    PHYSIK_API int PHYSIK_GetSurfaceTriangleIndexCount(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle surfaceExtraction);

    PHYSIK_API int PHYSIK_CopySurfaceTriangleIndices(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle surfaceExtraction,
        int* outIndices,
        int maxIndexCount);

    PHYSIK_API PhysiK::ComponentHandle PHYSIK_CreateSurfaceVisualComponent(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle surfaceExtractionHandle);

    PHYSIK_API int PHYSIK_GetSurfaceVisualVertexCount(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle surfaceVisualHandle);

    PHYSIK_API int PHYSIK_GetSurfaceVisualTriangleIndexCount(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle surfaceVisualHandle);

    PHYSIK_API int PHYSIK_GetSurfaceVisualNormalCount(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle surfaceVisualHandle);

    PHYSIK_API int PHYSIK_GetSurfaceVisualVertex(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle surfaceVisualHandle,
        int visualVertexIndex,
        float* outX,
        float* outY,
        float* outZ);

    PHYSIK_API int PHYSIK_GetSurfaceVisualTriangleIndex(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle surfaceVisualHandle,
        int triangleIndexArrayIndex,
        int* outIndex);

    PHYSIK_API int PHYSIK_GetSurfaceVisualNormal(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle surfaceVisualHandle,
        int visualNormalIndex,
        float* outX,
        float* outY,
        float* outZ);

    PHYSIK_API PhysiK::ComponentHandle PHYSIK_CreateTetMeshMapperComponent(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle sourceTetMesh,
        PhysiK::ComponentHandle destinationTetMesh);

    PHYSIK_API int PHYSIK_SetTetMeshLocalCurrentPosition(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle tetMesh,
        int localNodeIndex,
        float x,
        float y,
        float z);

    PHYSIK_API int PHYSIK_GetTetMeshLocalCurrentPosition(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle tetMesh,
        int localNodeIndex,
        float* outX,
        float* outY,
        float* outZ);

    PHYSIK_API PhysiK::ComponentHandle PHYSIK_CreateVisualMeshComponent(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle hostTetMesh);

    PHYSIK_API void PHYSIK_SetVisualMeshData(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle visualMesh,
        const PhysiK::Vec3* vertices,
        int vertexCount,
        const int* triangleIndices,
        int triangleIndexCount);

    PHYSIK_API int PHYSIK_BuildVisualMeshEmbedding(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle visualMesh);

    PHYSIK_API int PHYSIK_GetVisualMeshVertexCount(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle visualMesh);

    PHYSIK_API int PHYSIK_GetVisualMeshTriangleIndexCount(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle visualMesh);

    PHYSIK_API int PHYSIK_CopyVisualMeshVertices(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle visualMesh,
        PhysiK::Vec3* outVertices,
        int maxVertexCount);

    PHYSIK_API int PHYSIK_CopyVisualMeshTriangleIndices(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle visualMesh,
        int* outIndices,
        int maxIndexCount);

    PHYSIK_API void PHYSIK_SetCollisionSphereConnectionSettings(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle sphereComponent,
        float stiffness,
        float damping);

    PHYSIK_API void PHYSIK_GetCollisionSphereConnectionSettings(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle sphereComponent,
        float* outStiffness,
        float* outDamping);

    PHYSIK_API void PHYSIK_DestroyComponent(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle component);

    PHYSIK_API int PHYSIK_IsComponentHandleValid(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle component);

    PHYSIK_API int PHYSIK_GetTetMeshTetCount(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle component);

    PHYSIK_API int PHYSIK_GetTetMeshNodeCount(
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
