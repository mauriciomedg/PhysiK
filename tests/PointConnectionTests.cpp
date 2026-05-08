#include "PhysiK/API/PhysiKAPI.h"

#include <cassert>
#include <cmath>

namespace
{
    struct Point
    {
        float x = 0.0f;
        float y = 0.0f;
        float z = 0.0f;
    };

    Point GetNodePosition(PhysiK::WorldHandle world, int nodeIndex)
    {
        Point point;
        PHYSIK_GetNodePosition(world, nodeIndex, &point.x, &point.y, &point.z);
        return point;
    }

    Point BarycentricPoint(
        const Point& p0,
        const Point& p1,
        const Point& p2,
        const Point& p3,
        float w0,
        float w1,
        float w2,
        float w3)
    {
        return Point{
            p0.x * w0 + p1.x * w1 + p2.x * w2 + p3.x * w3,
            p0.y * w0 + p1.y * w1 + p2.y * w2 + p3.y * w3,
            p0.z * w0 + p1.z * w1 + p2.z * w2 + p3.z * w3};
    }

    float DistanceSquared(const Point& a, const Point& b)
    {
        const float dx = a.x - b.x;
        const float dy = a.y - b.y;
        const float dz = a.z - b.z;
        return dx * dx + dy * dy + dz * dz;
    }

    void CreateSingleTet(PhysiK::WorldHandle world, int (&outNodes)[4])
    {
        outNodes[0] = PHYSIK_AddNode(world, 0.0f, 0.0f, 0.0f, 1.0f);
        outNodes[1] = PHYSIK_AddNode(world, 1.0f, 0.0f, 0.0f, 1.0f);
        outNodes[2] = PHYSIK_AddNode(world, 0.0f, 1.0f, 0.0f, 1.0f);
        outNodes[3] = PHYSIK_AddNode(world, 0.0f, 0.0f, 1.0f, 1.0f);
        PHYSIK_AddTet(world, outNodes[0], outNodes[1], outNodes[2], outNodes[3]);

        const int tets[] = {0};
        assert(PHYSIK_CreateTetMeshComponent(world, outNodes, 4, tets, 1) != nullptr);
    }

    Point GetTetCentroid(PhysiK::WorldHandle world, const int (&nodes)[4])
    {
        return BarycentricPoint(
            GetNodePosition(world, nodes[0]),
            GetNodePosition(world, nodes[1]),
            GetNodePosition(world, nodes[2]),
            GetNodePosition(world, nodes[3]),
            0.25f,
            0.25f,
            0.25f,
            0.25f);
    }

    struct ExternalLogicTestState
    {
        int callbackCount = 0;
        int nodes[4] = {};
    };

    void AddPointConnectionFromExternalLogic(PhysiK::WorldHandle world, void* userData)
    {
        ExternalLogicTestState* state = static_cast<ExternalLogicTestState*>(userData);
        ++state->callbackCount;

        PHYSIK_AddPointConnection(
            world,
            state->nodes[0],
            state->nodes[1],
            state->nodes[2],
            state->nodes[3],
            0.25f,
            0.25f,
            0.25f,
            0.25f,
            0.25f,
            0.25f,
            1.25f,
            100.0f,
            0.0f);
    }
}

void ManualPointConnectionMovesBarycentricPoint()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    int nodes[4] = {};
    CreateSingleTet(world, nodes);

    const float w0 = 0.25f;
    const float w1 = 0.25f;
    const float w2 = 0.25f;
    const float w3 = 0.25f;
    const Point target{0.25f, 0.25f, 1.25f};

    const Point before = BarycentricPoint(
        GetNodePosition(world, nodes[0]),
        GetNodePosition(world, nodes[1]),
        GetNodePosition(world, nodes[2]),
        GetNodePosition(world, nodes[3]),
        w0,
        w1,
        w2,
        w3);

    PHYSIK_AddPointConnection(
        world,
        nodes[0],
        nodes[1],
        nodes[2],
        nodes[3],
        w0,
        w1,
        w2,
        w3,
        target.x,
        target.y,
        target.z,
        100.0f,
        0.0f);

    PHYSIK_Step(world, 0.1f);

    const Point after = BarycentricPoint(
        GetNodePosition(world, nodes[0]),
        GetNodePosition(world, nodes[1]),
        GetNodePosition(world, nodes[2]),
        GetNodePosition(world, nodes[3]),
        w0,
        w1,
        w2,
        w3);

    assert(DistanceSquared(after, target) < DistanceSquared(before, target));
    assert(after.z > before.z);

    PHYSIK_DestroyWorld(world);
}

void SphereContactCreatesTransientConnectionAndMovesTet()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    int nodes[4] = {};
    CreateSingleTet(world, nodes);

    assert(PHYSIK_CreateCollisionSphereComponent(world, 0.25f, 0.25f, 0.25f, 0.75f) != nullptr);
    assert(PHYSIK_GetPointConnectionCount(world) == 0);

    const Point before = GetTetCentroid(world, nodes);

    PHYSIK_Step(world, 0.1f);

    const Point after = GetTetCentroid(world, nodes);

    assert(after.z > before.z);
    assert(PHYSIK_GetPointConnectionCount(world) == 0);

    PHYSIK_DestroyWorld(world);
}

void ExternalLogicHookRunsOnceBeforeSubsteps()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    ExternalLogicTestState state;
    CreateSingleTet(world, state.nodes);
    PHYSIK_SetSubstepCount(world, 4);
    PHYSIK_SetExternalLogicCallback(world, AddPointConnectionFromExternalLogic, &state);

    const Point before = GetTetCentroid(world, state.nodes);

    PHYSIK_Step(world, 0.1f);

    const Point after = GetTetCentroid(world, state.nodes);

    assert(state.callbackCount == 1);
    assert(after.z > before.z);
    assert(PHYSIK_GetPointConnectionCount(world) == 0);

    PHYSIK_ClearExternalLogicCallback(world);
    PHYSIK_DestroyWorld(world);
}



void KinematicUpdateRunsAfterExternalLogicBeforePhysicsSubsteps()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    int nodes[4] = {};
    CreateSingleTet(world, nodes);

    const PhysiK::ComponentHandle sphere =
        PHYSIK_CreateCollisionSphereComponent(world, 2.0f, 2.0f, 2.0f, 0.75f);
    assert(sphere != nullptr);

    struct State
    {
        PhysiK::ComponentHandle sphere = nullptr;
        int callbackCount = 0;
    } state{sphere, 0};

    auto moveSphereInExternalLogic = [](PhysiK::WorldHandle worldHandle, void* userData)
    {
        auto* callbackState = static_cast<State*>(userData);
        ++callbackState->callbackCount;
        PHYSIK_SetCollisionComponentKinematicTarget(worldHandle, callbackState->sphere, 0.25f, 0.25f, 0.25f);
    };

    PHYSIK_SetExternalLogicCallback(world, moveSphereInExternalLogic, &state);
    PHYSIK_SetSubstepCount(world, 3);

    const Point before = GetTetCentroid(world, nodes);
    PHYSIK_Step(world, 0.1f);
    const Point after = GetTetCentroid(world, nodes);

    assert(state.callbackCount == 1);
    assert(after.z > before.z);

    PHYSIK_ClearExternalLogicCallback(world);
    PHYSIK_DestroyWorld(world);
}

int main()
{
    ManualPointConnectionMovesBarycentricPoint();
    SphereContactCreatesTransientConnectionAndMovesTet();
    ExternalLogicHookRunsOnceBeforeSubsteps();
    KinematicUpdateRunsAfterExternalLogicBeforePhysicsSubsteps();
    return 0;
}
