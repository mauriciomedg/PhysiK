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

- World
- Components
- Bodies
- Collision primitives
- Meshes
- Materials
- Physics models owned by components
- Component-local topology

Transient objects:

- Contacts
- Connections
- Solver energies
- Temporary solver data

Most important rule:

Only the solver modifies simulation state.

Collision systems and gameplay systems only generate contacts, events, or transient connections.

Connections are not components.

Connections are not persistent simulation objects.

Connections are transient penalty-energy terms generated for one frame or one physics substep.

---

# Ownership Philosophy

World owns global simulation state and orchestration.

Components own their private simulation structure and physics model.

World owns:

- global nodes
- component storage
- component handles
- collision detection engine
- solver systems
- transient connections
- simulation pipeline
- substepping

World does not own:

- tetrahedra
- triangles
- rod elements
- FEM element topology
- cloth topology
- rigid-body model internals
- component-owned physics models
- concrete component construction

Tetrahedra are an FEM abstraction.

Therefore, tetrahedra belong to TetMeshComponent or the FEM model owned by TetMeshComponent.

Triangles are a cloth/shell abstraction.

Therefore, triangles belong to TriMeshComponent or the ClothModel owned by TriMeshComponent.

Rod elements are a Cosserat rod abstraction.

Therefore, rod elements belong to LineMeshComponent or the CosseratRodModel owned by LineMeshComponent.

Rigid-body state belongs to RigidBodyComponent.

The World only stores global nodes so different systems can assemble into one solver system.

---

# Component-Owned Models

Components own the physics models that define their behavior.

Component creation logic belongs to component factories or the C API layer.

World registers already-created components.

World does not construct component internals.

Examples:

TetMeshComponent owns:

- tetrahedral topology
- material
- FEMModel

RigidBodyComponent owns:

- rigid-body state
- RigidBodyModel

LineMeshComponent owns:

- line or rod topology
- CosseratRodModel

TriMeshComponent owns:

- triangle topology
- ClothModel

CollisionComponent owns:

- collision shape data
- collision parameters

CollisionComponent talks to CollisionDetectionEngine.

CollisionComponent may generate contacts or transient physics connections.

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

- Own global nodes.
- Own all components.
- Own collision engine.
- Own solver systems.
- Run simulation pipeline.
- Manage substepping.
- Store transient connections.
- Clear transient connections after substeps.
- Export simulation state to Unity/Unreal/host application.

World does not own tetrahedra.

World does not own triangle topology.

World does not own rod topology.

World does not own component-local physics models.

World should not become a container for every possible physics topology.

class World
{
public:
    void Step(float frameDt);
    ComponentHandle AddComponent(
        std::unique_ptr<Component> component);

private:
    std::vector<Node> nodes;

    std::vector<std::unique_ptr<Component>> components;
    std::vector<std::unique_ptr<PhysicsConnection>> transientConnections;

    CollisionDetectionEngine collisionDetectionEngine;

    int substepCount = 1;
};

Important:

Do not add this to World:

std::vector<Tet> tets;

Tets are owned by TetMeshComponent.

Do not add this to World:

std::vector<PhysicsModel*> physicsModels;

Physics models are owned by components.

World calls Component::UpdateSystem, and each component delegates to its owned model.

World should orchestrate through Component and PhysicsConnection abstractions.

World should expose generic registration:

ComponentHandle AddComponent(std::unique_ptr<Component> component);

World node creation is geometric/topological:

int AddNode(const Vec3& position);

Use explicit fixed-state APIs for anchored nodes:

void SetNodeFixed(int nodeIndex, bool fixed);
bool IsNodeFixed(int nodeIndex) const;

Node inverse mass is assembled into SolverData for the current solve.

Public FEM workflows should not use inverse mass as a mass tuning API.

Do not add component-specific creation methods to World, such as:

World::CreateTetMeshComponent(...);

World::CreateCollisionSphereComponent(...);

Concrete component construction belongs to component factories/static Create methods or the C API layer.

Do not add typed component registries for orchestration, such as:

std::vector<TetMeshComponent*> tetMeshes;

std::vector<CollisionComponent*> collisionComponents;

Do not add typed transient connection arrays for orchestration, such as:

std::vector<PointConnection> pointConnections;

Transient connections are stored through PhysicsConnection pointers.

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
Component::UpdateFrame(...)
    ↓
Animation / kinematic target update
    - driven by gameplay state
    - updates animated targets
    - updates kinematic components
Component::UpdateKinematicTarget(...)
    ↓
For each physics substep:
        Collision detection
        Contact generation
        Convert collision contacts into transient PhysicsConnections
        SolverData clear
        Component::QueryContacts(...)
        Component::UpdateSystem(...)
            - TetMeshComponent calls FEMModel::UpdateSystem(...)
            - RigidBodyComponent calls RigidBodyModel::UpdateSystem(...)
            - LineMeshComponent calls CosseratRodModel::UpdateSystem(...)
            - TriMeshComponent calls ClothModel::UpdateSystem(...)
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

Components, component-owned physics models, and physics connections write their contributions into SolverData.

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

- Components do not solve physics.
- Component-owned physics models contribute internal material behavior.
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

Components own their private simulation data.

Components own their physics model when they have physical behavior.

Components expose UpdateSystem so World can ask them to contribute to SolverData.

## Base Component

class Component
{
public:
    bool active = true;

    virtual void UpdateFrame(World& world, float dt) {}

    virtual void UpdateKinematicTarget(World& world) {}

    virtual void QueryContacts(
        World& world,
        CollisionDetectionEngine& collisionDetectionEngine,
        std::vector<Contact>& outContacts) {}

    virtual void UpdateSystem(
        World& world,
        SolverData& solverData,
        float dt) {}

    virtual ~Component() = default;
};

Components do not solve physics.

Components do not directly modify simulation state.

Components delegate mathematical assembly to their owned physics model.

World should not use dynamic_cast to orchestrate component behavior.

Component behavior must be exposed through virtual hooks.

---

# TetMeshComponent

Volumetric deformable body.

TetMeshComponent owns the FEM topology.

TetMeshComponent owns its FEMModel.

TetMeshComponent owns the tetrahedra generated from a Unity/Unreal mesh.

When a mesh is tetrahedralized from Unity or Unreal:

- TetMeshComponent creates its local tetrahedral topology.
- TetMeshComponent adds nodes to World.
- TetMeshComponent stores global node indices.
- TetMeshComponent keeps tetrahedra private.
- World does not publish or own the tetrahedra.

class TetMeshComponent : public Component
{
public:
    static std::unique_ptr<TetMeshComponent> CreateFromGlobalNodes(
        World& world,
        const int* globalNodeIndices,
        int nodeCount,
        const int* tetGlobalNodeIndices,
        int tetCount,
        const Material& material);

    static std::unique_ptr<TetMeshComponent> CreateFromPositions(
        World& world,
        const Vec3* positions,
        const int* fixedNodeFlags,
        int nodeCount,
        const int* tetLocalNodeIndices,
        int tetCount,
        const Material& material);

    std::vector<int> nodeIndices;
    std::vector<Tet> tets;

    Material material;

    FEMModel femModel;

    void UpdateSystem(
        World& world,
        SolverData& solverData,
        float dt) override
    {
        femModel.UpdateSystem(
            world,
            *this,
            solverData,
            dt);
    }
};

Uses:

Core/Physics/FEM

Important:

Tets are not global World data.

Tets are private FEM topology.

TetMeshComponent factories create and initialize tetrahedra.

TetMeshComponent or FEMModel initializes FEM rest data.

World does not initialize tetrahedra or FEM rest data.

Other systems should not depend on World owning tets.

Other systems that need tet queries should go through TetMeshComponent or CollisionDetectionEngine.

---

# RigidBodyComponent

Persistent rigid-body simulation object.

RigidBodyComponent owns its rigid-body state.

RigidBodyComponent owns its RigidBodyModel.

class RigidBodyComponent : public Component
{
public:
    Vec3 position;
    Quaternion orientation;

    Vec3 linearVelocity;
    Vec3 angularVelocity;

    float mass;
    Mat3 inertiaTensor;

    RigidBodyModel rigidBodyModel;

    void UpdateSystem(
        World& world,
        SolverData& solverData,
        float dt) override
    {
        rigidBodyModel.UpdateSystem(
            world,
            *this,
            solverData,
            dt);
    }
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

LineMeshComponent owns its line/rod topology.

LineMeshComponent owns its CosseratRodModel.

class LineMeshComponent : public Component
{
public:
    std::vector<int> nodeIndices;

    std::vector<RodElement> rodElements;

    CosseratRodModel cosseratRodModel;

    void UpdateSystem(
        World& world,
        SolverData& solverData,
        float dt) override
    {
        cosseratRodModel.UpdateSystem(
            world,
            *this,
            solverData,
            dt);
    }
};

Uses:

Core/Physics/CosseratRod

---

# TriMeshComponent

Represents:

- cloth
- shell simulation
- thin surface simulation

TriMeshComponent owns its triangle topology.

TriMeshComponent owns its ClothModel.

class TriMeshComponent : public Component
{
public:
    std::vector<int> nodeIndices;
    std::vector<TriangleElement> triangles;

    ClothModel clothModel;

    void UpdateSystem(
        World& world,
        SolverData& solverData,
        float dt) override
    {
        clothModel.UpdateSystem(
            world,
            *this,
            solverData,
            dt);
    }
};

Uses:

Core/Physics/Cloth

---

# CollisionComponent

Persistent collision component.

Collision components:

- own collision shape data
- talk to the CollisionDetectionEngine
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

    void UpdateKinematicTarget(World& world) override;

    void QueryContacts(
        World& world,
        CollisionDetectionEngine& collisionDetectionEngine,
        std::vector<Contact>& outContacts) override;
};

CollisionComponent reads/talks to CollisionDetectionEngine.

CollisionDetectionEngine performs the broad phase and narrow phase work.

CollisionComponent decides whether detected contacts become events or transient physics connections.

Collision kinematic target application is handled through Component::UpdateKinematicTarget.

---

# Collision Components

class CollisionSphereComponent : public CollisionComponent
{
public:
    static std::unique_ptr<CollisionSphereComponent> Create(
        const Vec3& position,
        float radius);

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
    void UpdateFrame(World& world, float dt) override;
};

Examples:

- ManualGrabComponent
- MetricsLoggerComponent
- ToolInteractionComponent

---

# Collision Detection Engine

The collision engine only detects contacts.

It does not create solver responses automatically.

It does not directly deform simulation objects.

It can query component-owned topology through components.

class CollisionDetectionEngine
{
public:
    CollisionDetectionEngine() = default;

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

The broad phase may index:

- component bounds
- tetrahedra from TetMeshComponent
- surface triangles from TriMeshComponent
- collision primitives

But the broad phase does not own the topology.

It only references or caches acceleration data.

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

Can store references to:

- tetrahedra owned by TetMeshComponent
- surface triangles owned by TriMeshComponent
- collision primitives owned by CollisionComponent

SpatialHashMap does not own simulation topology.

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

The narrow phase can read TetMeshComponent private topology through controlled access/query methods.

It should not require World to own tetrahedra.

---

# Contacts

Contacts are information.

Contacts are not solver constraints automatically.

struct Contact
{
    ComponentHandle sourceComponent;
    ComponentHandle targetComponent;

    int node0;
    int node1;
    int node2;
    int node3;

    Vec4 barycentric;

    Vec3 worldPoint;
    Vec3 normal;

    float penetrationDepth;
};

The querying component decides the response.

A collision component can convert a contact into a transient connection.

For soft-body FEM contacts, the four node indices refer to global World node indices.

The tetrahedron that produced the contact is owned by TetMeshComponent, not World.

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

- the 4 global node indices of the tetrahedron
- barycentric coordinates defining the material point inside the tetrahedron
- target data
- stiffness
- damping

The connection does not need to own or access the Tet object.

The Tet object remains private to TetMeshComponent.

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

Gradient:

g_i = k * w_i * (p - targetPosition)

Force:

f_i = -g_i

Implicit stiffness:

K_ij = k * w_i * w_j * I3

PointConnection contributes both force and positive stiffness blocks into SolverData.

The implicit solver consumes those stiffness blocks with the same sign convention as FEM:

A += dtÂ² K

b += dt f - dtÂ² K v_current

Force and stiffness are distributed to the four tet nodes using barycentric weights.

Damping is currently assembled as a force contribution only; add a solver-side damping matrix before treating it implicitly.

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

Physics models are owned by components.

Physics models are not components.

Physics models are not globally owned by World.

Each physics model contributes to SolverData through UpdateSystem.

## Base PhysicsModel

class PhysicsModel
{
public:
    virtual ~PhysicsModel() = default;
};

The base class can remain minimal.

Each concrete model can expose an UpdateSystem overload that receives the owning component type.

Example:

class FEMModel : public PhysicsModel
{
public:
    void UpdateSystem(
        World& world,
        TetMeshComponent& owner,
        SolverData& solverData,
        float dt);
};

This allows FEMModel to access the TetMeshComponent that owns the tetrahedra.

---

# FEMModel

FEMModel is the physics model responsible for volumetric deformable simulation.

FEMModel is owned by TetMeshComponent.

FEMModel works on the tetrahedra owned by TetMeshComponent.

FEMModel does not ask World for a global tet list.

class FEMModel : public PhysicsModel
{
public:
    void UpdateSystem(
        World& world,
        TetMeshComponent& owner,
        SolverData& solverData,
        float dt);
};

Responsibilities:

- Iterate over owner.tets.
- Use owner.nodeIndices and global World nodes.
- Compute deformation gradients.
- Compute strain energy.
- Compute internal forces.
- Compute stiffness contributions.
- Write contributions into SolverData.

FEMModel does not directly move nodes.

FEMModel does not own World nodes.

FEMModel does not require World to own tets.

---

# Linear Tetrahedral FEM Assembly Status

The engine currently supports small-strain linear tetrahedral FEM assembly:

- epsilon = B u_e
- sigma = D epsilon
- f_int = V B^T sigma
- K_e = V B^T D B

Current limitations:

- Linear small-strain model only.
- Temporary per-node damping only.
- Explicit force integration path currently.
- Stiffness blocks are assembled for a future implicit solve.
- No corotational model yet.
- No Neo-Hookean model yet.
- No implicit Euler yet.
- No fracture or tearing yet.

Safety requirements:

- Degenerate tetrahedra must be skipped safely.
- Material parameters must not produce NaN forces or stiffness blocks.
- FEMModel must not move nodes directly.
- FEMModel must assemble through SolverData.
- Tets remain owned by TetMeshComponent.

---

# FEM Lumped Mass

TetMeshComponent computes simple lumped nodal mass after tet rest data is initialized.

For each tetrahedron:

tetMass = material.density * tet.restVolume

Each tet node receives:

nodeMass += tetMass / 4

Then, for non-fixed nodes, SolverData receives:

nodeMass += tetMass / 4

For shared nodes, mass contributions accumulate from connected tetrahedra.

Nodes marked fixed with SetNodeFixed remain fixed and receive no dynamic node mass.

For FEM tet mesh nodes, density-derived mass is the default physical mass.

Positive manually supplied inverse mass values are legacy convenience only, not FEM mass tuning.

Unity-facing FEM usage:

- Create nodes using PHYSIK_AddNode(world, x, y, z).
- Use PHYSIK_SetNodeFixed only for fixed nodes.
- Use material density for FEM mass.
- Author FEM material in Unity as a ScriptableObject.
- Transfer material across the DLL boundary with PhysikMaterialDesc.

World nodes do not own FEM lumped mass or inverse mass.

SolverData owns the assembled mass values consumed by integration and implicit solves.

---

# Material Authoring Boundary

Material is component-owned physical data.

Unity-facing code should not expose C++ classes or namespaces directly.

The DLL boundary uses a C-compatible descriptor:

struct PhysikMaterialDesc
{
    float density;
    float youngModulus;
    float poissonRatio;
    float damping;
};

The C API converts PhysikMaterialDesc into PhysiK::Material internally.

TetMeshComponent stores the internal Material and updates its private tetrahedra when material coefficients change.

World nodes must not store material, mass, inverse mass, or FEM tuning data.

---

# Implicit Euler FEM Solve

The engine has a first implicit Euler solve path for node-based FEM dynamics.

Target system:

A Δv = b

A = M + dt² K

b = dt f_total

After solving:

v_new = v_old + Δv

x_new = x_old + dt v_new

Current limitations:

- First implicit milestone uses a dense linear solve.
- No Newton iterations yet.
- No nonlinear FEM yet.
- No corotational FEM yet.
- No Neo-Hookean FEM yet.
- Damping matrix C is not fully implemented yet.
- Dynamic-fixed stiffness coupling is not moved to the right-hand side yet.
- Dense solver is acceptable for small tests.
- A future milestone should replace the dense solve with sparse CG or LDLT.

Architecture:

- FEMModel assembles forces and positive stiffness blocks.
- PhysicsConnections may assemble forces and stiffness blocks.
- SolverData stores force contributions and stiffness blocks.
- World/Solver builds the global node system.
- Only solver/integration code updates node positions and velocities.
- Components and physics models do not move nodes directly.

---

# RigidBodyModel

RigidBodyModel is the physics model responsible for rigid-body dynamics.

RigidBodyModel is owned by RigidBodyComponent.

class RigidBodyModel : public PhysicsModel
{
public:
    void UpdateSystem(
        World& world,
        RigidBodyComponent& owner,
        SolverData& solverData,
        float dt);
};

RigidBodyModel does not directly move rigid bodies.

RigidBodyModel contributes rigid-body force, mass, inertia, and Jacobian terms into SolverData.

---

# CosseratRodModel

CosseratRodModel is the physics model responsible for rods, sutures, cables, catheters, and beams.

CosseratRodModel is owned by LineMeshComponent.

class CosseratRodModel : public PhysicsModel
{
public:
    void UpdateSystem(
        World& world,
        LineMeshComponent& owner,
        SolverData& solverData,
        float dt);
};

CosseratRodModel does not directly move nodes.

CosseratRodModel reads rod topology from the owning LineMeshComponent.

---

# ClothModel

ClothModel is the physics model responsible for cloth, shells, and thin surface simulation.

ClothModel is owned by TriMeshComponent.

class ClothModel : public PhysicsModel
{
public:
    void UpdateSystem(
        World& world,
        TriMeshComponent& owner,
        SolverData& solverData,
        float dt);
};

ClothModel does not directly move nodes.

ClothModel reads triangle topology from the owning TriMeshComponent.

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
    = scene-facing persistent objects that own private topology and models

Physics models
    = mathematical behavior owned by components

Collision engine
    = contact detection and spatial queries

Contacts
    = information

PhysicsConnections
    = transient solver energies

Solvers
    = numerical algorithms

World
    = global node storage and orchestration

---

# Most Important Rules

Collision systems do not deform.

Gameplay systems do not deform.

Components do not directly deform.

Physics models do not directly deform.

Physics connections do not directly deform.

Only the solver modifies simulation state.

World does not own tetrahedra.

World does not own component-local topology.

World does not globally own physics models.

Components own their physics models.

TetMeshComponent owns tetrahedra.

FEMModel works on the TetMeshComponent that owns it.

---

# DLL API

Expose stable C ABI.

extern "C"
{
    struct PhysikMaterialDesc
    {
        float density;
        float youngModulus;
        float poissonRatio;
        float damping;
    };

    PHYSIK_API WorldHandle PHYSIK_CreateWorld();

    PHYSIK_API void PHYSIK_DestroyWorld(
        WorldHandle world);

    PHYSIK_API void PHYSIK_Step(
        WorldHandle world,
        float dt);

    PHYSIK_API int PHYSIK_AddNode(
        WorldHandle world,
        float x,
        float y,
        float z);

    PHYSIK_API void PHYSIK_SetNodeFixed(
        WorldHandle world,
        int nodeIndex,
        int fixed);

    PHYSIK_API int PHYSIK_IsNodeFixed(
        WorldHandle world,
        int nodeIndex);

    PHYSIK_API ComponentHandle PHYSIK_CreateTetMeshComponent(...);

    PHYSIK_API ComponentHandle PHYSIK_CreateTetMeshComponentWithMaterialDesc(
        WorldHandle world,
        const int* nodeIndices,
        int nodeCount,
        const int* tetNodeIndices,
        int tetCount,
        const PhysikMaterialDesc* material);

    PHYSIK_API void PHYSIK_SetTetMeshMaterial(
        WorldHandle world,
        ComponentHandle component,
        const PhysikMaterialDesc* material);

    PHYSIK_API ComponentHandle PHYSIK_CreateRigidBodyComponent(...);

    PHYSIK_API ComponentHandle PHYSIK_CreateLineMeshComponent(...);

    PHYSIK_API ComponentHandle PHYSIK_CreateTriMeshComponent(...);

    PHYSIK_API ComponentHandle PHYSIK_CreateCollisionSphereComponent(...);

    PHYSIK_API int PHYSIK_GetTetMeshTetCount(...);

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

PHYSIK_CreateTetMeshComponent receives global node indices and tetrahedron node-index tuples.
The C API calls TetMeshComponent::CreateFromGlobalNodes and registers the result with World::AddComponent.

The created TetMeshComponent owns the tetrahedra.

Do not route new code through a global PHYSIK_AddTet-style World tet store.

Do not add a World::AddTet function.

Do not add World::CreateTetMeshComponent or World::CreateCollisionSphereComponent.

The C API can know which concrete component to create; World should not.

---

# Example World Step

void World::Step(float frameDt)
{
    ExternalSyncFromHost();

    for (std::unique_ptr<Component>& component : components)
    {
        if (component && component->active)
        {
            component->UpdateFrame(*this, frameDt);
        }
    }

    UpdateAnimationTargets(frameDt);

    float subDt = frameDt / static_cast<float>(substepCount);

    for (int substep = 0; substep < substepCount; ++substep)
    {
        std::vector<Contact> contacts;

        for (std::unique_ptr<Component>& component : components)
        {
            if (component && component->active)
            {
                component->QueryContacts(
                    *this,
                    collisionDetectionEngine,
                    contacts);
            }
        }

        GenerateCollisionConnections();

        SolverData solverData;
        solverData.Clear();

        for (std::unique_ptr<Component>& component : components)
        {
            if (component && component->active)
            {
                component->UpdateSystem(
                    *this,
                    solverData,
                    subDt);
            }
        }

        for (std::unique_ptr<PhysicsConnection>& connection : transientConnections)
        {
            connection->UpdateSystem(*this, solverData, subDt);
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
- global nodes
- component storage
- SolverData
- Component::UpdateSystem
- Component::UpdateFrame
- Component::UpdateKinematicTarget
- Component::QueryContacts
- TetMeshComponent
- TetMeshComponent-owned tetrahedra
- TetMeshComponent-owned FEMModel placeholder
- PhysicsConnections folder
- PhysicsConnection base class
- PointConnection
- transientConnections stored as PhysicsConnection pointers
- simple solver
- substepping
- transient connection clearing
- DLL API

Goal:

- One TetMeshComponent
- One tetrahedron owned by TetMeshComponent
- Global nodes stored in World
- No global World tet array
- One transient PointConnection
- One solver step
- One substep clear

This validates the ownership architecture.

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

2. Respect the ownership rule: World owns global nodes and orchestration, components own topology and physics models.

3. Do not make World own tetrahedra.

4. Do not make World own triangle topology.

5. Do not make World own rod topology.

6. Do not make World globally own PhysicsModel instances.

7. TetMeshComponent owns its tetrahedra.

8. TetMeshComponent owns its FEMModel.

9. RigidBodyComponent owns its RigidBodyModel.

10. LineMeshComponent owns its CosseratRodModel.

11. TriMeshComponent owns its ClothModel.

12. World calls Component::UpdateSystem during solver assembly.

13. Each component delegates UpdateSystem to its owned model.

14. Do not turn PointConnection, SurfaceConnection, or LineConnection into components.

15. Do not make PointConnection, SurfaceConnection, or LineConnection persistent scene objects.

16. Store connections under Core/PhysicsConnections.

17. Connections must expose UpdateSystem.

18. Physics models must contribute to SolverData through their owning component.

19. Physics connections contribute to SolverData.

20. Collision code must not directly modify node positions or velocities.

21. Gameplay code must not directly modify node positions or velocities.

22. Components must not directly solve physics.

23. World orchestrates the simulation pipeline.

24. The solver is the only system allowed to modify simulation state.

25. World must use Component hooks instead of dynamic_cast-based orchestration.

26. Collision-generated contacts are converted into transient PhysicsConnections before connection assembly.

27. World::Solve is currently a temporary explicit-force solve path that transfers SolverData nodal forces into Node::force.

28. Keep the first milestone minimal.

29. Do not over-engineer future systems before the first vertical slice is working.

---

# Final Architecture Reminder

World says:

I own global nodes and orchestrate the simulation.

A TetMeshComponent says:

This deformable object exists in the scene, and I own its tetrahedral topology.

A FEMModel says:

Given my owning TetMeshComponent, I know how to assemble deformable FEM equations.

A PointConnection says:

For this substep, this material point defined by four global nodes and barycentric coordinates should be pulled toward this target.

A solver says:

Given all model energies and connection energies, compute the next physical state.

That is the core architecture of PhysiK.
