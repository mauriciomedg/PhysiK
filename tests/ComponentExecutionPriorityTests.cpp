#include "PhysiK/Components/Component.h"
#include "PhysiK/Components/ComponentExecutionPriority.h"
#include "PhysiK/Components/ScriptComponent.h"
#include "PhysiK/Components/SurfaceExtractionComponent.h"
#include "PhysiK/Components/SurfaceVisualComponent.h"
#include "PhysiK/Components/TetMeshComponent.h"
#include "PhysiK/Components/TetMeshMapperComponent.h"
#include "PhysiK/Components/TetMeshPhysicsComponent.h"
#include "PhysiK/Components/TopologyMeshComponent.h"
#include "PhysiK/Components/VisualMeshComponent.h"
#include "PhysiK/Core/PhysicsConnections/PhysicsConnection.h"
#include "PhysiK/Core/Solvers/SolverData.h"
#include "PhysiK/Core/World/World.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <memory>
#include <vector>

namespace
{
    class DefaultPriorityTestComponent :
        public PhysiK::Component
    {
    };

    class RecordingComponent :
        public PhysiK::Component
    {
    public:
        RecordingComponent(
            PhysiK::ComponentExecutionPriority priority,
            int id,
            std::vector<int>& callOrder)
            : priority(priority)
            , id(id)
            , callOrder(callOrder)
        {
        }

        PhysiK::ComponentExecutionPriority
        GetExecutionPriority() const override
        {
            return priority;
        }

        void PreUpdate(
            PhysiK::World&,
            float) override
        {
            callOrder.push_back(id);
        }

    private:
        PhysiK::ComponentExecutionPriority priority;
        int id;
        std::vector<int>& callOrder;
    };

    class RecordingConnection :
        public PhysiK::PhysicsConnection
    {
    public:
        explicit RecordingConnection(int& assemblyCount)
            : assemblyCount(assemblyCount)
        {
        }

        void UpdateSystem(
            PhysiK::World&,
            PhysiK::SolverData&,
            float) override
        {
            ++assemblyCount;
        }

    private:
        int& assemblyCount;
    };

    class FrameConnectionComponent :
        public PhysiK::Component
    {
    public:
        FrameConnectionComponent(
            int& preUpdateCount,
            int& frameConnectionAssemblyCount)
            : preUpdateCount(preUpdateCount)
            , frameConnectionAssemblyCount(frameConnectionAssemblyCount)
        {
        }

        void PreUpdate(
            PhysiK::World& world,
            float) override
        {
            ++preUpdateCount;
            world.AddTransientConnection(
                std::make_unique<RecordingConnection>(
                    frameConnectionAssemblyCount));
        }

    private:
        int& preUpdateCount;
        int& frameConnectionAssemblyCount;
    };

    class SubstepConnectionComponent :
        public PhysiK::Component
    {
    public:
        SubstepConnectionComponent(
            int& createdCount,
            int& substepConnectionAssemblyCount)
            : createdCount(createdCount)
            , substepConnectionAssemblyCount(substepConnectionAssemblyCount)
        {
        }

        void UpdateSystem(
            PhysiK::World& world,
            PhysiK::SolverData&,
            float) override
        {
            ++createdCount;
            world.AddTransientConnection(
                std::make_unique<RecordingConnection>(
                    substepConnectionAssemblyCount));
        }

    private:
        int& createdCount;
        int& substepConnectionAssemblyCount;
    };

    void Require(bool condition, const char* message)
    {
        if (!condition)
        {
            std::fprintf(stderr, "ComponentExecutionPriorityTests failed: %s\n", message);
            std::exit(1);
        }
    }

    std::unique_ptr<RecordingComponent> MakeRecordingComponent(
        PhysiK::ComponentExecutionPriority priority,
        int id,
        std::vector<int>& callOrder)
    {
        return std::make_unique<RecordingComponent>(
            priority,
            id,
            callOrder);
    }
}

void DefaultComponentPriorityIsDefault()
{
    DefaultPriorityTestComponent component;

    Require(
        component.GetExecutionPriority() ==
            PhysiK::ComponentExecutionPriority::
                Default,
        "default component priority should be Default");
}

void ComponentExecutionPriorityLessOrdersValues()
{
    const PhysiK::ComponentExecutionPriorityLess less;

    Require(
        less(
            PhysiK::ComponentExecutionPriority::
                ScriptComponent,
            PhysiK::ComponentExecutionPriority::
                CollisionSphereComponent),
        "ScriptComponent should sort before CollisionSphereComponent");
    Require(
        !less(
            PhysiK::ComponentExecutionPriority::
                VisualMeshComponent,
            PhysiK::ComponentExecutionPriority::
                ScriptComponent),
        "VisualMeshComponent should not sort before ScriptComponent");
}

void ConcreteComponentsReportExecutionPriorities()
{
    using PhysiK::ComponentExecutionPriority;

    PhysiK::ScriptComponent scriptComponent;
    PhysiK::TetMeshComponent tetMeshComponent;
    PhysiK::TetMeshPhysicsComponent tetMeshPhysicsComponent;
    PhysiK::TopologyMeshComponent topologyMeshComponent;
    PhysiK::TetMeshMapperComponent tetMeshMapperComponent;
    PhysiK::SurfaceExtractionComponent surfaceExtractionComponent;
    PhysiK::SurfaceVisualComponent surfaceVisualComponent;
    PhysiK::VisualMeshComponent visualMeshComponent;

    Require(
        scriptComponent.GetExecutionPriority() ==
            ComponentExecutionPriority::ScriptComponent,
        "ScriptComponent should report ScriptComponent priority");
    Require(
        tetMeshComponent.GetExecutionPriority() ==
            ComponentExecutionPriority::TetMeshComponent,
        "TetMeshComponent should report TetMeshComponent priority");
    Require(
        tetMeshPhysicsComponent.GetExecutionPriority() ==
            ComponentExecutionPriority::TetMeshPhysicsComponent,
        "TetMeshPhysicsComponent should report TetMeshPhysicsComponent priority");
    Require(
        topologyMeshComponent.GetExecutionPriority() ==
            ComponentExecutionPriority::TopologyMeshComponent,
        "TopologyMeshComponent should report TopologyMeshComponent priority");
    Require(
        tetMeshMapperComponent.GetExecutionPriority() ==
            ComponentExecutionPriority::TetMeshMapperComponent,
        "TetMeshMapperComponent should report TetMeshMapperComponent priority");
    Require(
        surfaceExtractionComponent.GetExecutionPriority() ==
            ComponentExecutionPriority::SurfaceExtractionComponent,
        "SurfaceExtractionComponent should report SurfaceExtractionComponent priority");
    Require(
        surfaceVisualComponent.GetExecutionPriority() ==
            ComponentExecutionPriority::SurfaceVisualComponent,
        "SurfaceVisualComponent should report SurfaceVisualComponent priority");
    Require(
        visualMeshComponent.GetExecutionPriority() ==
            ComponentExecutionPriority::VisualMeshComponent,
        "VisualMeshComponent should report VisualMeshComponent priority");
}

void ComponentExecutionPriorityLessOrdersAdjacentHierarchy()
{
    using PhysiK::ComponentExecutionPriority;

    const PhysiK::ComponentExecutionPriorityLess less;
    const std::vector<ComponentExecutionPriority> hierarchy{
        ComponentExecutionPriority::ScriptComponent,
        ComponentExecutionPriority::CollisionSphereComponent,
        ComponentExecutionPriority::TetMeshComponent,
        ComponentExecutionPriority::TetMeshPhysicsComponent,
        ComponentExecutionPriority::TopologyMeshComponent,
        ComponentExecutionPriority::TetMeshMapperComponent,
        ComponentExecutionPriority::SurfaceExtractionComponent,
        ComponentExecutionPriority::SurfaceVisualComponent,
        ComponentExecutionPriority::VisualMeshComponent,
        ComponentExecutionPriority::Default};

    for (std::size_t index = 1; index < hierarchy.size(); ++index)
    {
        Require(
            less(
                hierarchy[index - 1u],
                hierarchy[index]),
            "component execution priorities should sort adjacent hierarchy entries");
    }
}

void ComponentExecutionPriorityMultimapPreservesPriorityOrderAndDuplicates()
{
    using PhysiK::ComponentExecutionPriority;

    std::multimap<
        ComponentExecutionPriority,
        int,
        PhysiK::ComponentExecutionPriorityLess>
        orderedValues;

    orderedValues.insert({ComponentExecutionPriority::Default, 10});
    orderedValues.insert({ComponentExecutionPriority::TetMeshPhysicsComponent, 4});
    orderedValues.insert({ComponentExecutionPriority::CollisionSphereComponent, 2});
    orderedValues.insert({ComponentExecutionPriority::VisualMeshComponent, 9});
    orderedValues.insert({ComponentExecutionPriority::ScriptComponent, 1});
    orderedValues.insert({ComponentExecutionPriority::SurfaceVisualComponent, 8});
    orderedValues.insert({ComponentExecutionPriority::TetMeshComponent, 3});
    orderedValues.insert({ComponentExecutionPriority::TetMeshMapperComponent, 6});
    orderedValues.insert({ComponentExecutionPriority::TopologyMeshComponent, 5});
    orderedValues.insert({ComponentExecutionPriority::SurfaceExtractionComponent, 7});
    orderedValues.insert({ComponentExecutionPriority::CollisionSphereComponent, 20});

    std::vector<ComponentExecutionPriority> priorities;
    std::vector<int> collisionSphereValues;
    for (const auto& entry : orderedValues)
    {
        priorities.push_back(entry.first);
        if (entry.first == ComponentExecutionPriority::CollisionSphereComponent)
        {
            collisionSphereValues.push_back(entry.second);
        }
    }

    const std::vector<ComponentExecutionPriority> expectedPriorities{
        ComponentExecutionPriority::ScriptComponent,
        ComponentExecutionPriority::CollisionSphereComponent,
        ComponentExecutionPriority::CollisionSphereComponent,
        ComponentExecutionPriority::TetMeshComponent,
        ComponentExecutionPriority::TetMeshPhysicsComponent,
        ComponentExecutionPriority::TopologyMeshComponent,
        ComponentExecutionPriority::TetMeshMapperComponent,
        ComponentExecutionPriority::SurfaceExtractionComponent,
        ComponentExecutionPriority::SurfaceVisualComponent,
        ComponentExecutionPriority::VisualMeshComponent,
        ComponentExecutionPriority::Default};

    Require(
        priorities == expectedPriorities,
        "multimap iteration should follow component execution priority order");
    Require(
        collisionSphereValues.size() == 2u,
        "multimap should preserve duplicate CollisionSphereComponent entries");
    Require(
        collisionSphereValues[0] == 2 &&
            collisionSphereValues[1] == 20,
        "duplicate CollisionSphereComponent entries should preserve inserted values");
}

void WorldExecutesComponentsInPriorityOrder()
{
    using PhysiK::ComponentExecutionPriority;

    PhysiK::World world;
    std::vector<int> callOrder;

    world.AddComponent(
        MakeRecordingComponent(
            ComponentExecutionPriority::VisualMeshComponent,
            9,
            callOrder));
    world.AddComponent(
        MakeRecordingComponent(
            ComponentExecutionPriority::TetMeshPhysicsComponent,
            4,
            callOrder));
    world.AddComponent(
        MakeRecordingComponent(
            ComponentExecutionPriority::ScriptComponent,
            1,
            callOrder));
    world.AddComponent(
        MakeRecordingComponent(
            ComponentExecutionPriority::CollisionSphereComponent,
            2,
            callOrder));

    world.Step(0.0f);

    const std::vector<int> expectedOrder{1, 2, 4, 9};
    Require(
        callOrder == expectedOrder,
        "World should execute PreUpdate in component priority order");
}

void WorldExecutesDuplicatePriorityComponents()
{
    using PhysiK::ComponentExecutionPriority;

    PhysiK::World world;
    std::vector<int> callOrder;

    world.AddComponent(
        MakeRecordingComponent(
            ComponentExecutionPriority::CollisionSphereComponent,
            10,
            callOrder));
    world.AddComponent(
        MakeRecordingComponent(
            ComponentExecutionPriority::CollisionSphereComponent,
            20,
            callOrder));

    world.Step(0.0f);

    Require(
        callOrder.size() == 2u,
        "World should execute both duplicate-priority components");
    Require(
        (callOrder[0] == 10 && callOrder[1] == 20) ||
            (callOrder[0] == 20 && callOrder[1] == 10),
        "duplicate-priority components should preserve both recorded ids");
}

void WorldExecutesAllScriptPriorityComponentsBeforeCollision()
{
    using PhysiK::ComponentExecutionPriority;

    PhysiK::World world;
    std::vector<int> callOrder;

    world.AddComponent(
        MakeRecordingComponent(
            ComponentExecutionPriority::CollisionSphereComponent,
            20,
            callOrder));
    world.AddComponent(
        MakeRecordingComponent(
            ComponentExecutionPriority::ScriptComponent,
            1,
            callOrder));
    world.AddComponent(
        MakeRecordingComponent(
            ComponentExecutionPriority::ScriptComponent,
            2,
            callOrder));
    world.AddComponent(
        MakeRecordingComponent(
            ComponentExecutionPriority::ScriptComponent,
            3,
            callOrder));

    world.Step(0.0f);

    Require(
        callOrder.size() == 4u,
        "World should execute all ScriptComponent-priority and collision components");
    Require(
        callOrder[3] == 20,
        "CollisionSphereComponent-priority component should execute after all scripts");

    std::vector<int> scriptIds{
        callOrder[0],
        callOrder[1],
        callOrder[2]};
    std::sort(
        scriptIds.begin(),
        scriptIds.end());

    const std::vector<int> expectedScriptIds{1, 2, 3};
    Require(
        scriptIds == expectedScriptIds,
        "all ScriptComponent-priority instances should execute before collision");
}

void WorldUnregistersDestroyedComponentsFromExecution()
{
    using PhysiK::ComponentExecutionPriority;

    PhysiK::World world;
    std::vector<int> callOrder;

    world.AddComponent(
        MakeRecordingComponent(
            ComponentExecutionPriority::ScriptComponent,
            1,
            callOrder));
    const PhysiK::ComponentHandle destroyedHandle =
        world.AddComponent(
            MakeRecordingComponent(
                ComponentExecutionPriority::CollisionSphereComponent,
                2,
                callOrder));
    world.AddComponent(
        MakeRecordingComponent(
            ComponentExecutionPriority::VisualMeshComponent,
            3,
            callOrder));

    world.DestroyComponent(destroyedHandle);
    world.Step(0.0f);

    const std::vector<int> expectedOrder{1, 3};
    Require(
        callOrder == expectedOrder,
        "destroyed component should be absent from ordered execution");
}

void WorldReusedSlotRegistersNewComponentOnly()
{
    using PhysiK::ComponentExecutionPriority;

    PhysiK::World world;
    std::vector<int> callOrder;

    const PhysiK::ComponentHandle destroyedHandle =
        world.AddComponent(
            MakeRecordingComponent(
                ComponentExecutionPriority::CollisionSphereComponent,
                1,
                callOrder));

    world.DestroyComponent(destroyedHandle);

    world.AddComponent(
        MakeRecordingComponent(
            ComponentExecutionPriority::ScriptComponent,
            2,
            callOrder));

    world.Step(0.0f);

    const std::vector<int> expectedOrder{2};
    Require(
        callOrder == expectedOrder,
        "reused slot should execute only the newly registered component");
}

void WorldComponentHandlesRemainStableWithOrderedExecution()
{
    using PhysiK::ComponentExecutionPriority;

    PhysiK::World world;
    std::vector<int> callOrder;

    const PhysiK::ComponentHandle firstHandle =
        world.AddComponent(
            MakeRecordingComponent(
                ComponentExecutionPriority::VisualMeshComponent,
                1,
                callOrder));
    const PhysiK::ComponentHandle secondHandle =
        world.AddComponent(
            MakeRecordingComponent(
                ComponentExecutionPriority::ScriptComponent,
                2,
                callOrder));

    PhysiK::Component* firstComponent =
        world.GetComponent(firstHandle);
    PhysiK::Component* secondComponent =
        world.GetComponent(secondHandle);

    Require(
        firstComponent != nullptr,
        "first component handle should resolve");
    Require(
        secondComponent != nullptr,
        "second component handle should resolve");
    Require(
        firstComponent->GetExecutionPriority() ==
            ComponentExecutionPriority::VisualMeshComponent,
        "first component handle should still refer to the vector slot component");
    Require(
        secondComponent->GetExecutionPriority() ==
            ComponentExecutionPriority::ScriptComponent,
        "second component handle should still refer to the vector slot component");

    const PhysiK::ComponentHandle firstIndexHandle =
        world.GetComponentHandleByIndex(
            static_cast<int>(
                firstHandle.index));
    const PhysiK::ComponentHandle secondIndexHandle =
        world.GetComponentHandleByIndex(
            static_cast<int>(
                secondHandle.index));

    Require(
        firstIndexHandle.index == firstHandle.index &&
            firstIndexHandle.generation == firstHandle.generation,
        "GetComponentHandleByIndex should preserve first handle");
    Require(
        secondIndexHandle.index == secondHandle.index &&
            secondIndexHandle.generation == secondHandle.generation,
        "GetComponentHandleByIndex should preserve second handle");
}

void WorldStateAddsRotationalDefaults()
{
    PhysiK::WorldState state;

    const int stateIndex = state.AddNodeState(
        PhysiK::Vec3{1.0f, 2.0f, 3.0f},
        PhysiK::Vec3{4.0f, 5.0f, 6.0f},
        7.0f);

    Require(
        stateIndex == 0 && state.IsValidStateIndex(stateIndex),
        "WorldState should create a valid first state index");
    Require(
        state.positions.size() == 1u &&
            state.positions.size() == state.velocities.size() &&
            state.positions.size() == state.forces.size() &&
            state.positions.size() == state.masses.size() &&
            state.positions.size() == state.orientations.size() &&
            state.positions.size() == state.angularVelocities.size() &&
            state.positions.size() == state.torques.size() &&
            state.positions.size() == state.inverseInertias.size(),
        "WorldState arrays should stay size-synchronized");

    const PhysiK::Quaternion& orientation = state.orientations[0];
    Require(
        orientation.w == 1.0f &&
            orientation.x == 0.0f &&
            orientation.y == 0.0f &&
            orientation.z == 0.0f,
        "WorldState should default orientation to identity");
    Require(
        state.angularVelocities[0].x == 0.0f &&
            state.angularVelocities[0].y == 0.0f &&
            state.angularVelocities[0].z == 0.0f &&
            state.torques[0].x == 0.0f &&
            state.torques[0].y == 0.0f &&
            state.torques[0].z == 0.0f,
        "WorldState should default angular velocity and torque to zero");

    const PhysiK::Mat3& inverseInertia = state.inverseInertias[0];
    Require(
        inverseInertia.columns[0].x == 0.0f &&
            inverseInertia.columns[0].y == 0.0f &&
            inverseInertia.columns[0].z == 0.0f &&
            inverseInertia.columns[1].x == 0.0f &&
            inverseInertia.columns[1].y == 0.0f &&
            inverseInertia.columns[1].z == 0.0f &&
            inverseInertia.columns[2].x == 0.0f &&
            inverseInertia.columns[2].y == 0.0f &&
            inverseInertia.columns[2].z == 0.0f,
        "WorldState should default inverse inertia to zero");
}

void WorldNodesReferenceWorldStateStorage()
{
    PhysiK::World world;

    const int firstNode = world.AddNode(PhysiK::Vec3{1.0f, 2.0f, 3.0f});
    const int secondNode = world.AddNode(PhysiK::Vec3{4.0f, 5.0f, 6.0f});

    Require(
        firstNode == 0 && secondNode == 1,
        "node indices should remain sequential");
    Require(
        world.GetNode(firstNode).stateIndex == firstNode &&
            world.GetNode(secondNode).stateIndex == secondNode,
        "node metadata should point to matching world state slots");
    Require(
        !world.NodeHasRotation(firstNode) &&
            !world.NodeHasRotation(secondNode),
        "new FEM-style nodes should not advertise rotational state");

    Require(
        world.GetNodePosition(firstNode).x == 1.0f &&
            world.GetNodePosition(firstNode).y == 2.0f &&
            world.GetNodePosition(firstNode).z == 3.0f,
        "node position should be stored in WorldState");

    world.GetNodeVelocity(firstNode) = PhysiK::Vec3{7.0f, 8.0f, 9.0f};
    Require(
        world.GetNodeVelocity(firstNode).x == 7.0f &&
            world.GetNodeVelocity(firstNode).y == 8.0f &&
            world.GetNodeVelocity(firstNode).z == 9.0f,
        "node velocity accessor should reference WorldState");

    world.SetNodePosition(firstNode, PhysiK::Vec3{10.0f, 11.0f, 12.0f});
    Require(
        world.GetNodePosition(firstNode).x == 10.0f &&
            world.GetNodePosition(firstNode).y == 11.0f &&
            world.GetNodePosition(firstNode).z == 12.0f,
        "SetNodePosition should update WorldState position");
    Require(
        world.GetNodeVelocity(firstNode).x == 0.0f &&
            world.GetNodeVelocity(firstNode).y == 0.0f &&
            world.GetNodeVelocity(firstNode).z == 0.0f,
        "SetNodePosition should preserve existing velocity reset behavior");

    world.SetNodeHasRotation(firstNode, true);
    world.SetNodeOrientation(
        firstNode,
        PhysiK::Quaternion{0.5f, 0.5f, 0.5f, 0.5f});
    world.SetNodeAngularVelocity(firstNode, PhysiK::Vec3{1.0f, 2.0f, 3.0f});
    world.SetNodeTorque(firstNode, PhysiK::Vec3{4.0f, 5.0f, 6.0f});
    world.SetNodeInverseInertia(
        firstNode,
        PhysiK::Mat3::FromColumns(
            PhysiK::Vec3{1.0f, 0.0f, 0.0f},
            PhysiK::Vec3{0.0f, 2.0f, 0.0f},
            PhysiK::Vec3{0.0f, 0.0f, 3.0f}));

    Require(
        world.NodeHasRotation(firstNode),
        "rotational capability flag should be settable");
    Require(
        world.GetNodeOrientation(firstNode).w == 0.5f &&
            world.GetNodeOrientation(firstNode).x == 0.5f &&
            world.GetNodeOrientation(firstNode).y == 0.5f &&
            world.GetNodeOrientation(firstNode).z == 0.5f,
        "orientation accessor should reference WorldState");
    Require(
        world.GetNodeAngularVelocity(firstNode).x == 1.0f &&
            world.GetNodeAngularVelocity(firstNode).y == 2.0f &&
            world.GetNodeAngularVelocity(firstNode).z == 3.0f,
        "angular velocity accessor should reference WorldState");
    Require(
        world.GetNodeTorque(firstNode).x == 4.0f &&
            world.GetNodeTorque(firstNode).y == 5.0f &&
            world.GetNodeTorque(firstNode).z == 6.0f,
        "torque accessor should reference WorldState");
    Require(
        world.GetNodeInverseInertia(firstNode).columns[0].x == 1.0f &&
            world.GetNodeInverseInertia(firstNode).columns[1].y == 2.0f &&
            world.GetNodeInverseInertia(firstNode).columns[2].z == 3.0f,
        "inverse inertia accessor should reference WorldState");

    world.SetNodeFixed(firstNode, true);
    Require(
        world.GetNodeAngularVelocity(firstNode).x == 0.0f &&
            world.GetNodeAngularVelocity(firstNode).y == 0.0f &&
            world.GetNodeAngularVelocity(firstNode).z == 0.0f &&
            world.GetNodeTorque(firstNode).x == 0.0f &&
            world.GetNodeTorque(firstNode).y == 0.0f &&
            world.GetNodeTorque(firstNode).z == 0.0f,
        "fixed nodes should clear angular velocity and torque storage");
}

void WorldFrameConnectionsSurviveEverySubstep()
{
    PhysiK::World world;
    int preUpdateCount = 0;
    int frameConnectionAssemblies = 0;

    world.SetSubstepCount(3);
    world.AddComponent(
        std::make_unique<FrameConnectionComponent>(
            preUpdateCount,
            frameConnectionAssemblies));

    world.Step(0.1f);

    Require(
        preUpdateCount == 1,
        "frame connection component should run PreUpdate once per frame");
    Require(
        frameConnectionAssemblies == 3,
        "frame-scoped connection should assemble once in every substep");
    Require(
        world.GetTransientConnectionCount() == 0,
        "frame-scoped connections should be cleared after the frame");
}

void WorldSubstepConnectionsAreRecreatedAndTrimmed()
{
    PhysiK::World world;
    int substepConnectionsCreated = 0;
    int substepConnectionAssemblies = 0;

    world.SetSubstepCount(3);
    world.AddComponent(
        std::make_unique<SubstepConnectionComponent>(
            substepConnectionsCreated,
            substepConnectionAssemblies));

    world.Step(0.1f);

    Require(
        substepConnectionsCreated == 3,
        "substep component should create one connection per substep");
    Require(
        substepConnectionAssemblies == 3,
        "substep connections should assemble once and not survive into later substeps");
    Require(
        world.GetTransientConnectionCount() == 0,
        "substep-scoped connections should be cleared after the frame");
}

void WorldFrameAndSubstepConnectionsCoexist()
{
    PhysiK::World world;
    int preUpdateCount = 0;
    int frameConnectionAssemblies = 0;
    int substepConnectionsCreated = 0;
    int substepConnectionAssemblies = 0;

    world.SetSubstepCount(3);
    world.AddComponent(
        std::make_unique<FrameConnectionComponent>(
            preUpdateCount,
            frameConnectionAssemblies));
    world.AddComponent(
        std::make_unique<SubstepConnectionComponent>(
            substepConnectionsCreated,
            substepConnectionAssemblies));

    world.Step(0.1f);

    Require(
        preUpdateCount == 1,
        "frame connection component should run once when coexisting with substep connections");
    Require(
        frameConnectionAssemblies == 3,
        "frame connection should persist through all substeps");
    Require(
        substepConnectionsCreated == 3,
        "substep connection should be recreated each substep");
    Require(
        substepConnectionAssemblies == 3,
        "only the current substep connection should assemble each substep");
    Require(
        world.GetTransientConnectionCount() == 0,
        "all connections should be cleared after mixed frame/substep lifecycle");
}

void WorldClearsAllConnectionsAfterPositiveFrame()
{
    PhysiK::World world;
    int preUpdateCount = 0;
    int frameConnectionAssemblies = 0;

    world.SetSubstepCount(2);
    world.AddComponent(
        std::make_unique<FrameConnectionComponent>(
            preUpdateCount,
            frameConnectionAssemblies));

    world.Step(0.1f);

    Require(
        world.GetTransientConnectionCount() == 0,
        "positive-dt frame should flush the complete connection table");
}

void WorldZeroDtFrameFlushesFrameConnections()
{
    PhysiK::World world;
    int preUpdateCount = 0;
    int frameConnectionAssemblies = 0;

    world.AddComponent(
        std::make_unique<FrameConnectionComponent>(
            preUpdateCount,
            frameConnectionAssemblies));

    world.Step(0.0f);

    Require(
        preUpdateCount == 1,
        "zero-dt frame should still run PreUpdate");
    Require(
        frameConnectionAssemblies == 0,
        "zero-dt frame should not assemble frame connections");
    Require(
        world.GetTransientConnectionCount() == 0,
        "zero-dt frame should flush frame-scoped connection requests");
}

void WorldNegativeDtDoesNotRunConnectionLifecycle()
{
    PhysiK::World world;
    int preUpdateCount = 0;
    int frameConnectionAssemblies = 0;

    world.AddComponent(
        std::make_unique<FrameConnectionComponent>(
            preUpdateCount,
            frameConnectionAssemblies));

    world.Step(-1.0f);

    Require(
        preUpdateCount == 0,
        "negative dt should not run PreUpdate");
    Require(
        frameConnectionAssemblies == 0,
        "negative dt should not assemble connections");
    Require(
        world.GetTransientConnectionCount() == 0,
        "negative dt should not add or clear connections through lifecycle work");
}

void WorldOneSubstepConnectionLifecycleWorks()
{
    PhysiK::World world;
    int preUpdateCount = 0;
    int frameConnectionAssemblies = 0;
    int substepConnectionsCreated = 0;
    int substepConnectionAssemblies = 0;

    world.SetSubstepCount(1);
    world.AddComponent(
        std::make_unique<FrameConnectionComponent>(
            preUpdateCount,
            frameConnectionAssemblies));
    world.AddComponent(
        std::make_unique<SubstepConnectionComponent>(
            substepConnectionsCreated,
            substepConnectionAssemblies));

    world.Step(0.1f);

    Require(
        frameConnectionAssemblies == 1,
        "one-substep frame connection should assemble once");
    Require(
        substepConnectionsCreated == 1,
        "one-substep component should create one substep connection");
    Require(
        substepConnectionAssemblies == 1,
        "one-substep substep connection should assemble once");
    Require(
        world.GetTransientConnectionCount() == 0,
        "one-substep frame should clear the complete connection table");
}

int main()
{
    DefaultComponentPriorityIsDefault();
    ComponentExecutionPriorityLessOrdersValues();
    ConcreteComponentsReportExecutionPriorities();
    ComponentExecutionPriorityLessOrdersAdjacentHierarchy();
    ComponentExecutionPriorityMultimapPreservesPriorityOrderAndDuplicates();
    WorldExecutesComponentsInPriorityOrder();
    WorldExecutesDuplicatePriorityComponents();
    WorldExecutesAllScriptPriorityComponentsBeforeCollision();
    WorldUnregistersDestroyedComponentsFromExecution();
    WorldReusedSlotRegistersNewComponentOnly();
    WorldComponentHandlesRemainStableWithOrderedExecution();
    WorldStateAddsRotationalDefaults();
    WorldNodesReferenceWorldStateStorage();
    WorldFrameConnectionsSurviveEverySubstep();
    WorldSubstepConnectionsAreRecreatedAndTrimmed();
    WorldFrameAndSubstepConnectionsCoexist();
    WorldClearsAllConnectionsAfterPositiveFrame();
    WorldZeroDtFrameFlushesFrameConnections();
    WorldNegativeDtDoesNotRunConnectionLifecycle();
    WorldOneSubstepConnectionLifecycleWorks();

    return 0;
}
