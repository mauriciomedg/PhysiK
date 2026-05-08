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
- Physics models

Transient objects:

- Contacts
- Connections
- Solver energies

Most important rule:

Only the solver modifies simulation state.

Collision systems and gameplay systems only generate contacts, events, or transient connections.

Connections are not components.

Connections are not persistent simulation objects.

Connections are transient penalty-energy terms generated for one frame or one physics substep.

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
          PhysicsModel.h

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

        PhysicsConnections/
          PhysicsConnection.h

          PointConnection.h
          SurfaceConnection.h
          LineConnection.h

          RigidBodyConnection.h
          RigidBodyOrientationConnection.h

        Collision/
          CollisionDetectionEngine.h
          BroadPhase.h
          NarrowPhase.h
          SpatialHashMap.h

        Solvers/
          SolverData.h

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

  src/
    Core/
      World/
      Physics/
      PhysicsConnections/
      Collision/
      Solvers/

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
- Own all physics models.
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

    std::vector<PhysicsModel*> physicsModels;

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
        SolverData clear
        PhysicsModel::UpdateSystem(...)
        PhysicsConnection::UpdateSystem(...)
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

# SolverData

SolverData is the temporary structure used to assemble the current physics system.

Physics models and physics connections write their contributions into SolverData.

class SolverData
{
public:
    void Clear();

    // Force / residual vector
    // Stiffness matrix
    // Mass matrix
    // Damping matrix
    // Temporary solver buffers
};

Important:

- Physics models contribute internal material behavior.
- Physics connections contribute transient penalty energies.
- The solver consumes SolverData.
- Only the solver modifies the final simulation state.

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

Gameplay does not directly deform simulation objects.

Gameplay creates transient connections.

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

Components do not solve physics.

Components do not directly modify simulation state.

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
- ToolInteractionComponent

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

A collision component can convert a contact into a transient connection.

---

# Connection Architecture

Connections are transient solver inputs.

Connections are not components.

Connections are not persistent simulation objects.

Connections are consumed by the solver during one physics substep.

Connections disappear after the substep.

Only these systems create connections:

- CollisionComponent
- GameplayComponent
- internal physics interaction logic

Connections live in:

Core/PhysicsConnections/

A connection represents a penalty energy over simulation variables.

A connection contributes force, residual, stiffness, damping, or Jacobian terms into SolverData.

A connection does not directly modify:

- node positions
- node velocities
- rigid-body positions
- rigid-body orientations
- rigid-body velocities

## Base PhysicsConnection

class PhysicsConnection
{
public:
    virtual void UpdateSystem(
        World& world,
        SolverData& solverData,
        float dt) = 0;

    virtual ~PhysicsConnection() = default;
};

---

# Soft-Body Connections

Soft-body connections apply to a material point inside a tetrahedron.

The connection receives:

- the 4 nodes of the tetrahedron
- barycentric coordinates defining the material point inside the tetrahedron
- target data
- stiffness
- damping

Material point:

Vec3 p =
    barycentric.x * x0 +
    barycentric.y * x1 +
    barycentric.z * x2 +
    barycentric.w * x3;

---

## PointConnection

PointConnection pulls a tetrahedral material point toward a target world-space point.

struct PointConnection : public PhysicsConnection
{
    int node0;
    int node1;
    int node2;
    int node3;

    Vec4 barycentric;

    Vec3 targetPosition;

    float stiffness;
    float damping;

    void UpdateSystem(
        World& world,
        SolverData& solverData,
        float dt) override;
};

Energy:

E = 0.5 * k * ||p - targetPosition||²

Force is distributed to the four tet nodes using barycentric weights.

Use cases:

- grabbing
- tool interaction
- contact response
- manual gameplay connection
- temporary attachment

---

## SurfaceConnection

SurfaceConnection pushes or pulls a tetrahedral material point toward a plane or surface.

struct SurfaceConnection : public PhysicsConnection
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

    void UpdateSystem(
        World& world,
        SolverData& solverData,
        float dt) override;
};

Signed distance:

float d = dot(p - surfacePoint, surfaceNormal);

Energy:

E = 0.5 * k * d²

Use cases:

- plane contact
- wall contact
- surface projection
- tool-surface interaction

---

## LineConnection

LineConnection pulls a tetrahedral material point toward a line.

struct LineConnection : public PhysicsConnection
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

    void UpdateSystem(
        World& world,
        SolverData& solverData,
        float dt) override;
};

Closest point:

Vec3 closestPoint = ClosestPointOnLine(
    p,
    linePoint,
    lineDirection);

Energy:

E = 0.5 * k * ||p - closestPoint||²

Use cases:

- sliding interaction
- guide constraints
- catheter/tool alignment
- line-based temporary control

---

# Rigid Body Connections

No rigid-body constraints/joints yet.

Rigid-body interactions use transient energy connections.

## RigidBodyConnection

struct RigidBodyConnection : public PhysicsConnection
{
    RigidBodyHandle rigidBody;

    Vec3 localPoint;
    Vec3 targetPosition;

    float stiffness;
    float damping;

    void UpdateSystem(
        World& world,
        SolverData& solverData,
        float dt) override;
};

World-space point:

Vec3 p = rigidBody.position + rigidBody.orientation * localPoint;

Energy:

E = 0.5 * k * ||p - targetPosition||²

Produces:

- linear force
- angular torque

---

## RigidBodyOrientationConnection

struct RigidBodyOrientationConnection : public PhysicsConnection
{
    RigidBodyHandle rigidBody;

    Quaternion targetOrientation;

    float stiffness;
    float damping;

    void UpdateSystem(
        World& world,
        SolverData& solverData,
        float dt) override;
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
- stiffness matrices
- damping matrices
- mass matrices

Physics models do not hardcode numerical solvers.

Examples:

- Core/Physics/FEM
- Core/Physics/RigidBody
- Core/Physics/CosseratRod
- Core/Physics/Cloth

Physics models are persistent engine systems.

Physics models are not components.

Each physics model contributes to SolverData through UpdateSystem.

## Base PhysicsModel

class PhysicsModel
{
public:
    virtual void UpdateSystem(
        World& world,
        SolverData& solverData,
        float dt) = 0;

    virtual ~PhysicsModel() = default;
};

---

# FEMModel

FEMModel is the physics model responsible for volumetric deformable simulation.

class FEMModel : public PhysicsModel
{
public:
    void UpdateSystem(
        World& world,
        SolverData& solverData,
        float dt) override;
};

Responsibilities:

- Iterate over tetrahedra.
- Compute deformation gradients.
- Compute strain energy.
- Compute internal forces.
- Compute stiffness contributions.
- Write contributions into SolverData.

FEMModel does not directly move nodes.

---

# RigidBodyModel

RigidBodyModel is the physics model responsible for rigid-body dynamics.

class RigidBodyModel : public PhysicsModel
{
public:
    void UpdateSystem(
        World& world,
        SolverData& solverData,
        float dt) override;
};

RigidBodyModel does not directly move rigid bodies.

---

# CosseratRodModel

CosseratRodModel is the physics model responsible for rods, sutures, cables, catheters, and beams.

class CosseratRodModel : public PhysicsModel
{
public:
    void UpdateSystem(
        World& world,
        SolverData& solverData,
        float dt) override;
};

CosseratRodModel does not directly move nodes.

---

# ClothModel

ClothModel is the physics model responsible for cloth, shells, and thin surface simulation.

class ClothModel : public PhysicsModel
{
public:
    void UpdateSystem(
        World& world,
        SolverData& solverData,
        float dt) override;
};

ClothModel does not directly move nodes.

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
- Projective Dynamics

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

PhysicsConnections
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

Physics connections do not deform directly.

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

    PHYSIK_API void PHYSIK_AddPointConnection(
        WorldHandle world,
        int node0,
        int node1,
        int node2,
        int node3,
        float baryX,
        float baryY,
        float baryZ,
        float baryW,
        float targetX,
        float targetY,
        float targetZ,
        float stiffness,
        float damping);

    PHYSIK_API void PHYSIK_SetSubstepCount(
        WorldHandle world,
        int substepCount);
}

---

# Example World Step

void World::Step(float frameDt)
{
    ExternalSyncFromHost();

    for (GameplayComponent* gameplay : gameplayComponents)
    {
        if (gameplay && gameplay->active)
        {
            gameplay->Update(*this, frameDt);
        }
    }

    UpdateAnimationTargets(frameDt);

    float subDt = frameDt / static_cast<float>(substepCount);

    for (int substep = 0; substep < substepCount; ++substep)
    {
        std::vector<Contact> contacts;

        for (CollisionComponent* collision : collisionComponents)
        {
            if (collision && collision->active)
            {
                collision->QueryContacts(*this, contacts);
            }
        }

        GenerateConnectionsFromContacts(contacts);

        SolverData solverData;
        solverData.Clear();

        for (PhysicsModel* model : physicsModels)
        {
            if (model)
            {
                model->UpdateSystem(*this, solverData, subDt);
            }
        }

        for (PointConnection& connection : pointConnections)
        {
            connection.UpdateSystem(*this, solverData, subDt);
        }

        for (SurfaceConnection& connection : surfaceConnections)
        {
            connection.UpdateSystem(*this, solverData, subDt);
        }

        for (LineConnection& connection : lineConnections)
        {
            connection.UpdateSystem(*this, solverData, subDt);
        }

        for (RigidBodyConnection& connection : rigidBodyConnections)
        {
            connection.UpdateSystem(*this, solverData, subDt);
        }

        for (RigidBodyOrientationConnection& connection : rigidBodyOrientationConnections)
        {
            connection.UpdateSystem(*this, solverData, subDt);
        }

        solver.Solve(solverData, subDt);

        Integrate(subDt);

        ClearTransientConnections();
    }

    PostFrameUpdate();

    ExportDataToHost();
}

---

# Initial Milestone

First vertical slice:

- World
- SolverData
- PhysicsModel base class
- FEMModel placeholder
- TetMeshComponent
- PhysicsConnections folder
- PhysicsConnection base class
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

- full FEM
- collision detection
- rigid bodies
- rods
- cloth
- CG
- Unreal integration
- Unity integration
- tearing
- fracture
- topology changes

---

# Agent Instructions

When modifying this codebase:

1. Respect the separation between persistent simulation objects and transient solver inputs.

2. Do not turn PointConnection, SurfaceConnection, or LineConnection into components.

3. Do not make PointConnection, SurfaceConnection, or LineConnection persistent scene objects.

4. Store connections under Core/PhysicsConnections.

5. Connections must expose UpdateSystem.

6. Physics models must expose UpdateSystem.

7. Physics models and connections contribute to SolverData.

8. Collision code must not directly modify node positions or velocities.

9. Gameplay code must not directly modify node positions or velocities.

10. Components must not directly solve physics.

11. World orchestrates the simulation pipeline.

12. The solver is the only system allowed to modify simulation state.

13. Keep the first milestone minimal.

14. Do not over-engineer future systems before the first vertical slice is working.

---

# Final Architecture Reminder

A TetMeshComponent says:

This object exists in the scene.

A FEMModel says:

This object behaves like deformable FEM material.

A PointConnection says:

For this substep, this material point should be pulled toward this target.

A solver says:

Given all model energies and connection energies, compute the next physical state.

That is the core architecture of PhysiK.