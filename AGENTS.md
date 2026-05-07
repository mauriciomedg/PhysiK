# PhysiK Architecture

PhysiK is a modular C++ real-time physics simulation engine compiled as:

PhysiK.dll

Main namespace:

namespace PhysiK

DLL macro:

PHYSIK_API

---

# Core Philosophy

PhysiK separates:

- Persistent simulation objects
- Transient solver interactions

Persistent objects:

- Components
- Bodies
- Collision primitives
- Meshes
- Materials

Transient objects:

- Contacts
- Connections
- Solver energies

Most important rule:

Only the solver modifies simulation state.

Collision systems and gameplay systems only generate contacts, events, or transient connections.

---

# Folder Structure

PhysiK/
  include/
    PhysiK/

      API/
        PhysiKAPI.h
        Handles.h

      Core/
        World/
          World.h

        Physics/
          FEM/
            FEMModel.h
            TetElement.h
            FEMEnergy.h

          RigidBody/
            RigidBodyModel.h

          CosseratRod/
            CosseratRodModel.h
            RodElement.h
            RodEnergy.h

          Cloth/
            ClothModel.h
            TriangleElement.h
            ClothEnergy.h

        Collision/
          CollisionDetectionEngine.h
          BroadPhase.h
          NarrowPhase.h
          SpatialHashMap.h

        Solvers/
          Linear/
            ConjugateGradientSolver.h
            LDLTSolver.h
            GaussSeidelSolver.h

          Nonlinear/
            NewtonSolver.h
            ProjectiveDynamicsSolver.h

          Explicit/
            EulerIntegrator.h
            SemiImplicitEulerIntegrator.h
            TLEDSolver.h

      Components/
        Component.h

        TetMeshComponent.h
        RigidBodyComponent.h
        LineMeshComponent.h
        TriMeshComponent.h

        CollisionComponent.h
        CollisionSphereComponent.h
        CollisionCapsuleComponent.h
        CollisionBoxComponent.h
        CollisionPlaneComponent.h

        GameplayComponent.h

      Math/
        Vec3.h
        Vec4.h
        Mat3.h
        Quaternion.h
        Transform.h
        AABB.h

      PhysicsData/
        Node.h
        Tet.h
        Material.h

        Contact.h
        ContactCandidate.h

        PointConnection.h
        SurfaceConnection.h
        LineConnection.h

        RigidBodyConnection.h
        RigidBodyOrientationConnection.h

  src/
    Core/
    Components/
    Math/
    PhysicsData/

  tests/

---

# World

World is the main simulation container and orchestrator.

Responsibilities:

- Own simulation data.
- Own all components.
- Own collision engine.
- Own solver systems.
- Run simulation pipeline.
- Manage substepping.
- Store transient connections.
- Clear transient connections after substeps.

class World
{
public:
    void Step(float frameDt);

private:
    std::vector<Node> nodes;
    std::vector<Tet> tets;

    std::vector<Component*> components;

    std::vector<TetMeshComponent*> tetMeshes;
    std::vector<RigidBodyComponent*> rigidBodies;
    std::vector<LineMeshComponent*> lineMeshes;
    std::vector<TriMeshComponent*> triMeshes;

    std::vector<CollisionComponent*> collisionComponents;
    std::vector<GameplayComponent*> gameplayComponents;

    std::vector<PointConnection> pointConnections;
    std::vector<SurfaceConnection> surfaceConnections;
    std::vector<LineConnection> lineConnections;

    std::vector<RigidBodyConnection> rigidBodyConnections;
    std::vector<RigidBodyOrientationConnection> rigidBodyOrientationConnections;

    CollisionDetectionEngine collisionDetectionEngine;

    int substepCount = 1;
};

---

# Simulation Pipeline

Gameplay runs once per frame, before animation and before physics substeps.

Physics substeps are only for numerical stability.

World::Step(frameDt)
    ↓
External sync from Unity/Unreal
    ↓
GameplayComponent update
    - high-level logic
    - scripting logic
    - metrics
    - interaction commands
    - manual connection requests
    - simulation decisions
    ↓
Animation / kinematic target update
    - driven by gameplay state
    - updates animated targets
    - updates kinematic components
    ↓
For each physics substep:
        Collision detection
        Contact generation
        Internal connection generation
        Pre-update physics models
        Solver assembly
        Numerical solve
        Integration
        Post-substep physics update
        Clear transient substep connections
    ↓
Post-frame update
    ↓
Export updated data to Unity/Unreal

Important:

Connections live for one physics substep only.

If a contact persists, the connection must be recreated in the next substep.

Gameplay does not run inside the substep loop.

---

# Gameplay as Simulation Brain

The physics engine by itself is only a set of components, algorithms, models, and solvers.

GameplayComponent is the high-level brain that drives simulation behavior.

It decides what should happen during the frame:

- create interaction commands
- request manual connections
- enable/disable components
- drive animation targets
- trigger metrics/logging
- react to events
- define scripted physical behavior

Gameplay runs once per frame, after external sync and before animation/physics.

---

# Components

Components are persistent simulation objects attached to scene objects.

Components are exposed through the DLL.

## Base Component

class Component
{
public:
    bool active = true;

    virtual void Update(World& world, float dt) {}
    virtual ~Component() = default;
};

---

# TetMeshComponent

Volumetric deformable body.

class TetMeshComponent : public Component
{
public:
    std::vector<int> nodeIndices;
    std::vector<int> tetIndices;

    Material material;
};

Uses:

Core/Physics/FEM

---

# RigidBodyComponent

Persistent rigid-body simulation object.

class RigidBodyComponent : public Component
{
public:
    Vec3 position;
    Quaternion orientation;

    Vec3 linearVelocity;
    Vec3 angularVelocity;

    float mass;
    Mat3 inertiaTensor;
};

Uses:

Core/Physics/RigidBody

No persistent rigid-body constraints or joints for now.

Interactions are handled through transient rigid-body connections.

---

# LineMeshComponent

Represents:

- Cosserat rods
- beams
- sutures
- cables
- catheters

class LineMeshComponent : public Component
{
public:
    std::vector<int> nodeIndices;
};

Uses:

Core/Physics/CosseratRod

---

# TriMeshComponent

Represents:

- cloth
- shell simulation
- thin surface simulation

class TriMeshComponent : public Component
{
public:
    std::vector<int> nodeIndices;
    std::vector<int> triangleIndices;
};

Uses:

Core/Physics/Cloth

---

# CollisionComponent

Persistent collision component.

Collision components:

- query collision engine
- generate contacts
- optionally generate transient connections

They do not directly deform simulation objects.

class CollisionComponent : public Component
{
public:
    Transform transform;

    bool generateConnections = true;
    bool generateEvents = false;
    bool isSensor = false;

    float contactStiffness = 1000.0f;
    float contactDamping = 10.0f;

    virtual void QueryContacts(
        World& world,
        std::vector<Contact>& outContacts) = 0;
};

---

# Collision Components

class CollisionSphereComponent : public CollisionComponent
{
public:
    float radius;
};

class CollisionCapsuleComponent : public CollisionComponent
{
public:
    float radius;
    float height;
};

class CollisionBoxComponent : public CollisionComponent
{
public:
    Vec3 halfExtents;
};

class CollisionPlaneComponent : public CollisionComponent
{
public:
    Vec3 normal;
    float offset;
};

---

# GameplayComponent

Persistent scripting/high-level logic component.

Responsibilities:

- manual connection generation
- gameplay logic
- metrics
- logging
- event handling

class GameplayComponent : public Component
{
public:
    virtual void Update(World& world, float dt) = 0;
};

Examples:

- ManualGrabComponent
- MetricsLoggerComponent

---

# Collision Detection Engine

The collision engine only detects contacts.

It does not create solver responses automatically.

class CollisionDetectionEngine
{
public:
    void QueryContacts(
        World& world,
        const CollisionComponent& component,
        std::vector<Contact>& outContacts);

private:
    BroadPhase broadPhase;
    NarrowPhase narrowPhase;
};

---

# Broad Phase

Acceleration system.

Uses spatial hash map.

class BroadPhase
{
private:
    SpatialHashMap spatialHashMap;
};

Purpose:

Avoid N² collision checks.

---

# Spatial Hash Map

class SpatialHashMap
{
public:
    void Insert(...);
    void Query(...);

private:
    std::unordered_map<HashKey, Cell> cells;
};

Can store:

- tetrahedra
- surface triangles
- collision primitives

---

# Narrow Phase

Precise collision tests.

Examples:

- sphere vs tet
- sphere vs cloth
- capsule vs soft body
- box vs cloth
- plane vs FEM

Outputs:

Contact list

---

# Contacts

Contacts are information.

Contacts are not solver constraints automatically.

struct Contact
{
    int tetNode0;
    int tetNode1;
    int tetNode2;
    int tetNode3;

    Vec4 barycentric;

    Vec3 worldPoint;
    Vec3 normal;

    float penetrationDepth;
};

The querying component decides the response.

---

# Connection Architecture

Connections are transient solver inputs.

Connections are not components.

Connections are consumed by the solver during one physics substep.

Connections disappear after the substep.

Only these systems create connections:

- CollisionComponent
- GameplayComponent

---

# Soft-Body Connections

## PointConnection

struct PointConnection
{
    int node0;
    int node1;
    int node2;
    int node3;

    Vec4 barycentric;

    Vec3 targetPosition;

    float stiffness;
    float damping;
};

Material point:

Vec3 p =
    barycentric.x * x0 +
    barycentric.y * x1 +
    barycentric.z * x2 +
    barycentric.w * x3;

Energy:

E = 0.5 * k * ||p - target||²

Force is distributed to the four tet nodes using barycentric weights.

---

## SurfaceConnection

struct SurfaceConnection
{
    int node0;
    int node1;
    int node2;
    int node3;

    Vec4 barycentric;

    Vec3 surfacePoint;
    Vec3 surfaceNormal;

    float stiffness;
    float damping;
};

Energy:

d = dot(p - surfacePoint, surfaceNormal)
E = 0.5 * k * d²

---

## LineConnection

struct LineConnection
{
    int node0;
    int node1;
    int node2;
    int node3;

    Vec4 barycentric;

    Vec3 linePoint;
    Vec3 lineDirection;

    float stiffness;
    float damping;
};

Energy:

closestPoint = ClosestPointOnLine(p, linePoint, lineDirection)
E = 0.5 * k * ||p - closestPoint||²

---

# Rigid Body Connections

No rigid-body constraints/joints yet.

Rigid-body interactions use transient energy connections.

## RigidBodyConnection

struct RigidBodyConnection
{
    RigidBodyHandle rigidBody;

    Vec3 localPoint;
    Vec3 targetPosition;

    float stiffness;
    float damping;
};

World-space point:

Vec3 p = rigidBody.position + rigidBody.orientation * localPoint;

Energy:

E = 0.5 * k * ||p - target||²

Produces:

- linear force
- angular torque

---

## RigidBodyOrientationConnection

struct RigidBodyOrientationConnection
{
    RigidBodyHandle rigidBody;

    Quaternion targetOrientation;

    float stiffness;
    float damping;
};

Orientation error:

Quaternion qError = targetOrientation * inverse(currentOrientation);
Vec3 rotationError = Log(qError);

Energy:

E = 0.5 * k * dot(rotationError, rotationError)

Produces corrective torque.

---

# Physics Models

Physics models define:

- energies
- forces
- residuals
- Jacobians
- mass matrices

Physics models do not hardcode numerical solvers.

Examples:

- Core/Physics/FEM
- Core/Physics/RigidBody
- Core/Physics/CosseratRod
- Core/Physics/Cloth

---

# Solver Library

Solvers are reusable numerical algorithms.

Examples:

- CG
- LDLT
- Gauss-Seidel
- Newton
- TLED
- Euler
- Semi-Implicit Euler

Different physics models may use different solvers.

Do not hard-code one solver per physics model.

---

# Key Separation

Components
    = scene-facing persistent objects

Physics models
    = mathematical behavior

Collision engine
    = contact detection

Contacts
    = information

Connections
    = transient solver energies

Solvers
    = numerical algorithms

World
    = orchestration

---

# Most Important Rule

Collision systems do not deform.

Gameplay systems do not deform.

Components do not deform.

Only the solver modifies simulation state.

---

# DLL API

Expose stable C ABI.

extern "C"
{
    PHYSIK_API WorldHandle PHYSIK_CreateWorld();

    PHYSIK_API void PHYSIK_DestroyWorld(
        WorldHandle world);

    PHYSIK_API void PHYSIK_Step(
        WorldHandle world,
        float dt);

    PHYSIK_API ComponentHandle PHYSIK_CreateTetMeshComponent(...);

    PHYSIK_API ComponentHandle PHYSIK_CreateRigidBodyComponent(...);

    PHYSIK_API ComponentHandle PHYSIK_CreateLineMeshComponent(...);

    PHYSIK_API ComponentHandle PHYSIK_CreateTriMeshComponent(...);

    PHYSIK_API ComponentHandle PHYSIK_CreateCollisionSphereComponent(...);

    PHYSIK_API void PHYSIK_AddPointConnection(...);

    PHYSIK_API void PHYSIK_SetSubstepCount(...);
}

---

# Initial Milestone

First vertical slice:

- World
- TetMeshComponent
- PointConnection
- simple solver
- substepping
- transient connection clearing
- DLL API

Goal:

- One tetrahedron
- One transient PointConnection
- One solver step
- One substep clear

This validates the architecture.

Do not implement yet:

- FEM
- collision detection
- rigid bodies
- rods
- cloth
- CG
- Unreal integration
- Unity integration