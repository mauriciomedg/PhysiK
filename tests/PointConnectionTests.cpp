#include "PhysiK/API/PhysiKAPI.h"
#include "PhysiK/Core/Physics/FEM/FEMModel.h"
#include "PhysiK/Core/Solvers/SolverData.h"

#include <cassert>
#include <cmath>
#include <vector>

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

    Point GetNodeVelocity(PhysiK::WorldHandle world, int nodeIndex)
    {
        Point point;
        PHYSIK_GetNodeVelocity(world, nodeIndex, &point.x, &point.y, &point.z);
        return point;
    }

    float GetNodeInverseMass(PhysiK::WorldHandle world, int nodeIndex)
    {
        return PHYSIK_GetNodeInverseMass(world, nodeIndex);
    }

    bool IsNodeFixed(PhysiK::WorldHandle world, int nodeIndex)
    {
        return PHYSIK_IsNodeFixed(world, nodeIndex) != 0;
    }

    int AddNode(PhysiK::WorldHandle world, float x, float y, float z)
    {
        return PHYSIK_AddNode(world, x, y, z);
    }

    int AddFixedNode(PhysiK::WorldHandle world, float x, float y, float z)
    {
        const int node = PHYSIK_AddNode(world, x, y, z);
        PHYSIK_SetNodeFixed(world, node, 1);
        return node;
    }

    void SetNodeVelocity(PhysiK::WorldHandle world, int nodeIndex, const Point& velocity)
    {
        PHYSIK_SetNodeVelocity(world, nodeIndex, velocity.x, velocity.y, velocity.z);
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

    float LengthSquared(const PhysiK::Vec3& value)
    {
        return value.x * value.x + value.y * value.y + value.z * value.z;
    }

    bool IsFinite(float value)
    {
        return std::isfinite(value);
    }

    bool IsFinite(const PhysiK::Vec3& value)
    {
        return IsFinite(value.x) && IsFinite(value.y) && IsFinite(value.z);
    }

    bool IsFinite(const PhysiK::Mat3& matrix)
    {
        return IsFinite(matrix.columns[0]) &&
            IsFinite(matrix.columns[1]) &&
            IsFinite(matrix.columns[2]);
    }

    bool NearlyEqual(float a, float b, float tolerance = 0.00001f)
    {
        return std::abs(a - b) <= tolerance;
    }

    bool NearlyEqual(
        const PhysiK::Vec3& a,
        const PhysiK::Vec3& b,
        float tolerance = 0.00001f)
    {
        return NearlyEqual(a.x, b.x, tolerance) &&
            NearlyEqual(a.y, b.y, tolerance) &&
            NearlyEqual(a.z, b.z, tolerance);
    }

    float GetMat3Value(const PhysiK::Mat3& matrix, int row, int column)
    {
        const PhysiK::Vec3& sourceColumn = matrix.columns[column];
        if (row == 0)
        {
            return sourceColumn.x;
        }

        if (row == 1)
        {
            return sourceColumn.y;
        }

        return sourceColumn.z;
    }

    PhysiK::Vec3 SumForcesForNode(const PhysiK::SolverData& solverData, int node)
    {
        PhysiK::Vec3 total;
        for (const PhysiK::SolverData::NodeForce& force : solverData.GetNodeForces())
        {
            if (force.node == node)
            {
                total += force.force;
            }
        }
        return total;
    }

    const PhysiK::SolverData::StiffnessBlock* FindBlock(
        const PhysiK::SolverData& solverData,
        int nodeA,
        int nodeB)
    {
        for (const PhysiK::SolverData::StiffnessBlock& block : solverData.GetStiffnessBlocks())
        {
            if (block.nodeA == nodeA && block.nodeB == nodeB)
            {
                return &block;
            }
        }

        return nullptr;
    }

    std::vector<PhysiK::Node> CreateUnitTetNodes()
    {
        std::vector<PhysiK::Node> nodes(4);
        nodes[0].position = PhysiK::Vec3{0.0f, 0.0f, 0.0f};
        nodes[1].position = PhysiK::Vec3{1.0f, 0.0f, 0.0f};
        nodes[2].position = PhysiK::Vec3{0.0f, 1.0f, 0.0f};
        nodes[3].position = PhysiK::Vec3{0.0f, 0.0f, 1.0f};
        return nodes;
    }

    PhysiK::Tet CreateUnitTet(float youngModulus = 100.0f, float poissonRatio = 0.25f)
    {
        PhysiK::Tet tet;
        tet.node0 = 0;
        tet.node1 = 1;
        tet.node2 = 2;
        tet.node3 = 3;
        tet.youngModulus = youngModulus;
        tet.poissonRatio = poissonRatio;
        tet.damping = 0.0f;
        return tet;
    }

    void CreateSingleTet(PhysiK::WorldHandle world, int (&outNodes)[4])
    {
        outNodes[0] = AddNode(world, 0.0f, 0.0f, 0.0f);
        outNodes[1] = AddNode(world, 1.0f, 0.0f, 0.0f);
        outNodes[2] = AddNode(world, 0.0f, 1.0f, 0.0f);
        outNodes[3] = AddNode(world, 0.0f, 0.0f, 1.0f);

        const int tetNodeIndices[] = {outNodes[0], outNodes[1], outNodes[2], outNodes[3]};
        const PhysiK::ComponentHandle tetMesh =
            PHYSIK_CreateTetMeshComponentWithMaterial(
                world,
                outNodes,
                4,
                tetNodeIndices,
                1,
                24.0f,
                25.0f,
                0.3f,
                0.25f);
        assert(PHYSIK_IsComponentHandleValid(world, tetMesh) == 1);
    }

    void CreateSingleTetWithMaterial(
        PhysiK::WorldHandle world,
        int (&outNodes)[4],
        float youngModulus,
        float damping = 0.0f,
        float density = 24.0f)
    {
        outNodes[0] = AddNode(world, 0.0f, 0.0f, 0.0f);
        outNodes[1] = AddNode(world, 1.0f, 0.0f, 0.0f);
        outNodes[2] = AddNode(world, 0.0f, 1.0f, 0.0f);
        outNodes[3] = AddNode(world, 0.0f, 0.0f, 1.0f);

        const int tetNodeIndices[] = {outNodes[0], outNodes[1], outNodes[2], outNodes[3]};
        const PhysiK::ComponentHandle tetMesh =
            PHYSIK_CreateTetMeshComponentWithMaterial(
                world,
                outNodes,
                4,
                tetNodeIndices,
                1,
                density,
                youngModulus,
                0.3f,
                damping);
        assert(PHYSIK_IsComponentHandleValid(world, tetMesh) == 1);
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
    assert(PHYSIK_GetPointConnectionCount(world) == 1);

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

void GravityMovesDynamicNode()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    const int node = AddNode(world, 0.0f, 0.0f, 0.0f);
    PHYSIK_SetGravity(world, 0.0f, -10.0f, 0.0f);

    PHYSIK_Step(world, 0.1f);

    const Point after = GetNodePosition(world, node);
    assert(after.y < -0.09f);
    assert(after.y > -0.11f);

    PHYSIK_DestroyWorld(world);
}

void ImplicitEulerGravityMatchesSemiImplicitEulerForFreeNode()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    const int node = AddNode(world, 0.0f, 0.0f, 0.0f);
    PHYSIK_SetGravity(world, 0.0f, -10.0f, 0.0f);
    PHYSIK_SetSolverMode(world, 1);

    PHYSIK_Step(world, 0.1f);

    const Point after = GetNodePosition(world, node);
    assert(after.y < -0.09f);
    assert(after.y > -0.11f);

    PHYSIK_DestroyWorld(world);
}

void ImplicitEulerFixedNodeDoesNotMove()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    const int node = AddFixedNode(world, 0.0f, 1.0f, 0.0f);
    PHYSIK_SetGravity(world, 0.0f, -10.0f, 0.0f);
    PHYSIK_SetSolverMode(world, 1);

    const Point beforePosition = GetNodePosition(world, node);
    const Point beforeVelocity = GetNodeVelocity(world, node);

    PHYSIK_Step(world, 0.1f);

    const Point afterPosition = GetNodePosition(world, node);
    const Point afterVelocity = GetNodeVelocity(world, node);
    assert(DistanceSquared(afterPosition, beforePosition) < 0.000001f);
    assert(DistanceSquared(afterVelocity, beforeVelocity) < 0.000001f);

    PHYSIK_DestroyWorld(world);
}

void AddNodeCreatesDynamicGeometryNode()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    const int node = PHYSIK_AddNode(world, 0.0f, 0.0f, 0.0f);

    assert(!IsNodeFixed(world, node));
    assert(GetNodeInverseMass(world, node) > 0.0f);

    PHYSIK_DestroyWorld(world);
}

void SetNodeFixedControlsInternalInverseMass()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    const int node = PHYSIK_AddNode(world, 0.0f, 0.0f, 0.0f);

    PHYSIK_SetNodeFixed(world, node, 1);
    assert(IsNodeFixed(world, node));
    assert(GetNodeInverseMass(world, node) == 0.0f);

    PHYSIK_SetNodeFixed(world, node, 0);
    assert(!IsNodeFixed(world, node));
    assert(NearlyEqual(GetNodeInverseMass(world, node), 1.0f));

    PHYSIK_DestroyWorld(world);
}

void LegacyAddNodeWithInverseMassStillWorks()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    const int node = PHYSIK_AddNodeWithInverseMass(world, 0.0f, 0.0f, 0.0f, 0.5f);

    assert(!IsNodeFixed(world, node));
    assert(NearlyEqual(GetNodeInverseMass(world, node), 0.5f));

    PHYSIK_DestroyWorld(world);
}

void SphereContactCreatesTransientConnectionAndMovesTet()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    int nodes[4] = {};
    CreateSingleTet(world, nodes);

    const PhysiK::ComponentHandle sphere =
        PHYSIK_CreateCollisionSphereComponent(world, 0.25f, 0.25f, 0.25f, 0.75f);
    assert(PHYSIK_IsComponentHandleValid(world, sphere) == 1);
    assert(PHYSIK_GetPointConnectionCount(world) == 0);

    const Point before = GetTetCentroid(world, nodes);

    PHYSIK_Step(world, 0.1f);

    const Point after = GetTetCentroid(world, nodes);

    assert(after.z > before.z);
    assert(PHYSIK_GetPointConnectionCount(world) == 0);

    PHYSIK_DestroyWorld(world);
}

void MultipleForceSourcesCoexist()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    int nodes[4] = {};
    CreateSingleTet(world, nodes);

    PHYSIK_SetGravity(world, 0.0f, 0.0f, -1.0f);
    PHYSIK_AddPointConnection(
        world,
        nodes[0],
        nodes[1],
        nodes[2],
        nodes[3],
        0.25f,
        0.25f,
        0.25f,
        0.25f,
        0.25f,
        0.25f,
        1.25f,
        100.0f,
        0.0f);

    const Point before = GetTetCentroid(world, nodes);
    PHYSIK_Step(world, 0.1f);
    const Point after = GetTetCentroid(world, nodes);

    assert(after.z > before.z);
    assert(after.z < 1.25f);
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
    assert(PHYSIK_IsComponentHandleValid(world, sphere) == 1);

    struct State
    {
        PhysiK::ComponentHandle sphere;
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

void FEMElasticityMovesDistortedTetTowardRestShape()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    const int node0 = AddFixedNode(world, 0.0f, 0.0f, 0.0f);
    const int node1 = AddFixedNode(world, 1.0f, 0.0f, 0.0f);
    const int node2 = AddFixedNode(world, 0.0f, 1.0f, 0.0f);
    const int node3 = AddNode(world, 0.0f, 0.0f, 1.0f);

    const int nodes[] = {node0, node1, node2, node3};
    const int tetNodeIndices[] = {node0, node1, node2, node3};
    const PhysiK::ComponentHandle tetMesh =
        PHYSIK_CreateTetMeshComponentWithMaterial(
            world,
            nodes,
            4,
            tetNodeIndices,
            1,
            24.0f,
            25.0f,
            0.3f,
            0.25f);
    assert(PHYSIK_IsComponentHandleValid(world, tetMesh) == 1);

    const Point restPosition = GetNodePosition(world, node3);
    PHYSIK_SetNodePosition(world, node3, 0.0f, 0.0f, 1.25f);
    const Point distortedPosition = GetNodePosition(world, node3);

    PHYSIK_SetSubstepCount(world, 4);
    PHYSIK_Step(world, 0.1f);

    const Point after = GetNodePosition(world, node3);

    assert(DistanceSquared(after, restPosition) < DistanceSquared(distortedPosition, restPosition));
    assert(after.z < distortedPosition.z);

    PHYSIK_DestroyWorld(world);
}

void ImplicitEulerFEMTetMovesDistortedNodeTowardRestShape()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    const int node0 = AddFixedNode(world, 0.0f, 0.0f, 0.0f);
    const int node1 = AddFixedNode(world, 1.0f, 0.0f, 0.0f);
    const int node2 = AddFixedNode(world, 0.0f, 1.0f, 0.0f);
    const int node3 = AddNode(world, 0.0f, 0.0f, 1.0f);

    const int nodes[] = {node0, node1, node2, node3};
    const int tetNodeIndices[] = {node0, node1, node2, node3};
    const PhysiK::ComponentHandle tetMesh =
        PHYSIK_CreateTetMeshComponentWithMaterial(
            world,
            nodes,
            4,
            tetNodeIndices,
            1,
            24.0f,
            25.0f,
            0.3f,
            0.25f);
    assert(PHYSIK_IsComponentHandleValid(world, tetMesh) == 1);

    const Point restPosition = GetNodePosition(world, node3);
    PHYSIK_SetNodePosition(world, node3, 0.0f, 0.0f, 1.25f);
    const Point distortedPosition = GetNodePosition(world, node3);

    PHYSIK_SetSolverMode(world, 1);
    PHYSIK_Step(world, 0.25f);

    const Point after = GetNodePosition(world, node3);

    assert(IsFinite(after.x));
    assert(IsFinite(after.y));
    assert(IsFinite(after.z));
    assert(DistanceSquared(after, restPosition) < DistanceSquared(distortedPosition, restPosition));
    assert(after.z < distortedPosition.z);

    PHYSIK_DestroyWorld(world);
}

void ImplicitEulerAllDynamicTetResistsRestShapeVelocity()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    int nodes[4] = {};
    CreateSingleTetWithMaterial(world, nodes, 1000.0f);
    PHYSIK_SetSolverMode(world, 1);
    PHYSIK_SetGravity(world, 0.0f, 0.0f, 0.0f);
    SetNodeVelocity(world, nodes[3], Point{0.0f, 0.0f, 1.0f});

    PHYSIK_Step(world, 0.1f);

    const Point afterPosition = GetNodePosition(world, nodes[3]);
    const Point afterVelocity = GetNodeVelocity(world, nodes[3]);

    assert(IsFinite(afterPosition.z));
    assert(IsFinite(afterVelocity.z));
    assert(afterVelocity.z < 1.0f);
    assert(afterPosition.z < 1.1f);

    PHYSIK_DestroyWorld(world);
}

void ImplicitEulerPointAnchoredTetRecoversFreeNodes()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    int nodes[4] = {};
    CreateSingleTetWithMaterial(world, nodes, 1000.0f);
    PHYSIK_SetSolverMode(world, 1);
    PHYSIK_SetGravity(world, 0.0f, 0.0f, 0.0f);
    SetNodeVelocity(world, nodes[2], Point{0.0f, 1.0f, 0.0f});
    SetNodeVelocity(world, nodes[3], Point{0.0f, 0.0f, 1.0f});

    PHYSIK_AddPointConnection(
        world,
        nodes[0],
        nodes[0],
        nodes[0],
        nodes[0],
        1.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        10000.0f,
        0.0f);
    PHYSIK_AddPointConnection(
        world,
        nodes[1],
        nodes[1],
        nodes[1],
        nodes[1],
        1.0f,
        0.0f,
        0.0f,
        0.0f,
        1.0f,
        0.0f,
        0.0f,
        10000.0f,
        0.0f);

    PHYSIK_Step(world, 0.1f);

    const Point node2Velocity = GetNodeVelocity(world, nodes[2]);
    const Point node3Velocity = GetNodeVelocity(world, nodes[3]);

    assert(IsFinite(node2Velocity.y));
    assert(IsFinite(node3Velocity.z));
    assert(node2Velocity.y < 1.0f);
    assert(node3Velocity.z < 1.0f);

    PHYSIK_DestroyWorld(world);
}

void ImplicitEulerFEMRecoveryBeatsGravityOnlyMotion()
{
    PhysiK::WorldHandle gravityOnlyWorld = PHYSIK_CreateWorld();
    assert(gravityOnlyWorld != nullptr);
    const int gravityOnlyNode = AddNode(gravityOnlyWorld, 0.0f, 0.0f, 1.0f);
    PHYSIK_SetSolverMode(gravityOnlyWorld, 1);
    PHYSIK_SetGravity(gravityOnlyWorld, 0.0f, 0.0f, -1.0f);
    SetNodeVelocity(gravityOnlyWorld, gravityOnlyNode, Point{0.0f, 0.0f, 1.0f});
    PHYSIK_Step(gravityOnlyWorld, 0.1f);
    const Point gravityOnlyPosition = GetNodePosition(gravityOnlyWorld, gravityOnlyNode);
    PHYSIK_DestroyWorld(gravityOnlyWorld);

    PhysiK::WorldHandle femWorld = PHYSIK_CreateWorld();
    assert(femWorld != nullptr);
    int nodes[4] = {};
    CreateSingleTetWithMaterial(femWorld, nodes, 1000.0f);
    PHYSIK_SetSolverMode(femWorld, 1);
    PHYSIK_SetGravity(femWorld, 0.0f, 0.0f, -1.0f);
    SetNodeVelocity(femWorld, nodes[3], Point{0.0f, 0.0f, 1.0f});
    PHYSIK_Step(femWorld, 0.1f);
    const Point femPosition = GetNodePosition(femWorld, nodes[3]);
    const Point femVelocity = GetNodeVelocity(femWorld, nodes[3]);

    assert(IsFinite(femPosition.z));
    assert(IsFinite(femVelocity.z));
    assert(femPosition.z < gravityOnlyPosition.z);
    assert(femVelocity.z < 0.9f);

    PHYSIK_DestroyWorld(femWorld);
}

void FEMLumpedMassUsesDensityAndRestVolume()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    int nodes[4] = {};
    CreateSingleTetWithMaterial(world, nodes, 0.0f, 0.0f, 1.0f);

    const float expectedInverseMass = 24.0f;
    for (int node : nodes)
    {
        assert(NearlyEqual(GetNodeInverseMass(world, node), expectedInverseMass));
    }

    PHYSIK_DestroyWorld(world);
}

void FEMLumpedMassDensityReducesAccelerationFromSameForce()
{
    auto stepWithDensity = [](float density)
    {
        PhysiK::WorldHandle world = PHYSIK_CreateWorld();
        assert(world != nullptr);

        int nodes[4] = {};
        CreateSingleTetWithMaterial(world, nodes, 0.0f, 0.0f, density);

        PHYSIK_AddPointConnection(
            world,
            nodes[3],
            nodes[3],
            nodes[3],
            nodes[3],
            1.0f,
            0.0f,
            0.0f,
            0.0f,
            0.0f,
            0.0f,
            2.0f,
            1.0f,
            0.0f);
        PHYSIK_Step(world, 0.1f);

        const float velocity = GetNodeVelocity(world, nodes[3]).z;
        PHYSIK_DestroyWorld(world);
        return velocity;
    };

    const float lowDensityVelocity = stepWithDensity(1.0f);
    const float highDensityVelocity = stepWithDensity(4.0f);

    assert(lowDensityVelocity > highDensityVelocity);
    assert(highDensityVelocity > 0.0f);
}

void FEMLumpedMassPreservesFixedNodes()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    const int node0 = AddFixedNode(world, 0.0f, 0.0f, 0.0f);
    const int node1 = AddNode(world, 1.0f, 0.0f, 0.0f);
    const int node2 = AddNode(world, 0.0f, 1.0f, 0.0f);
    const int node3 = AddNode(world, 0.0f, 0.0f, 1.0f);
    const int nodes[] = {node0, node1, node2, node3};
    const int tetNodeIndices[] = {node0, node1, node2, node3};

    const PhysiK::ComponentHandle tetMesh =
        PHYSIK_CreateTetMeshComponentWithMaterial(
            world,
            nodes,
            4,
            tetNodeIndices,
            1,
            1.0f,
            0.0f,
            0.3f,
            0.0f);
    assert(PHYSIK_IsComponentHandleValid(world, tetMesh) == 1);

    assert(GetNodeInverseMass(world, node0) == 0.0f);
    assert(NearlyEqual(GetNodeInverseMass(world, node1), 24.0f));
    assert(NearlyEqual(GetNodeInverseMass(world, node2), 24.0f));
    assert(NearlyEqual(GetNodeInverseMass(world, node3), 24.0f));

    PHYSIK_DestroyWorld(world);
}

void FEMGravityAccelerationIsIndependentOfDensity()
{
    auto stepWithDensity = [](float density)
    {
        PhysiK::WorldHandle world = PHYSIK_CreateWorld();
        assert(world != nullptr);

        int nodes[4] = {};
        CreateSingleTetWithMaterial(world, nodes, 0.0f, 0.0f, density);
        PHYSIK_SetGravity(world, 0.0f, -10.0f, 0.0f);
        PHYSIK_Step(world, 0.1f);

        const Point position = GetNodePosition(world, nodes[3]);
        const Point velocity = GetNodeVelocity(world, nodes[3]);
        PHYSIK_DestroyWorld(world);
        return std::vector<float>{position.y, velocity.y};
    };

    const std::vector<float> lowDensityState = stepWithDensity(1.0f);
    const std::vector<float> highDensityState = stepWithDensity(4.0f);

    assert(NearlyEqual(lowDensityState[0], highDensityState[0], 0.0001f));
    assert(NearlyEqual(lowDensityState[1], highDensityState[1], 0.0001f));
    assert(lowDensityState[0] < 0.0f);
    assert(lowDensityState[1] < 0.0f);
}

void ImplicitEulerUsesStiffnessBlocks()
{
    PhysiK::WorldHandle noStiffnessWorld = PHYSIK_CreateWorld();
    assert(noStiffnessWorld != nullptr);
    PHYSIK_SetSolverMode(noStiffnessWorld, 1);
    const int freeNode = AddNode(noStiffnessWorld, 0.0f, 0.0f, 1.0f);

    PHYSIK_AddPointConnection(
        noStiffnessWorld,
        freeNode,
        freeNode,
        freeNode,
        freeNode,
        1.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        2.0f,
        10.0f,
        0.0f);
    PHYSIK_Step(noStiffnessWorld, 0.5f);
    const float velocityWithoutStiffness = GetNodeVelocity(noStiffnessWorld, freeNode).z;
    PHYSIK_DestroyWorld(noStiffnessWorld);

    PhysiK::WorldHandle stiffnessWorld = PHYSIK_CreateWorld();
    assert(stiffnessWorld != nullptr);
    PHYSIK_SetSolverMode(stiffnessWorld, 1);
    const int node0 = AddFixedNode(stiffnessWorld, 0.0f, 0.0f, 0.0f);
    const int node1 = AddFixedNode(stiffnessWorld, 1.0f, 0.0f, 0.0f);
    const int node2 = AddFixedNode(stiffnessWorld, 0.0f, 1.0f, 0.0f);
    const int node3 = AddNode(stiffnessWorld, 0.0f, 0.0f, 1.0f);
    const int componentNodes[] = {node0, node1, node2, node3};
    const int tetNodeIndices[] = {node0, node1, node2, node3};
    const PhysiK::ComponentHandle tetMesh =
        PHYSIK_CreateTetMeshComponentWithMaterial(
            stiffnessWorld,
            componentNodes,
            4,
            tetNodeIndices,
            1,
            24.0f,
            100.0f,
            0.3f,
            0.0f);
    assert(PHYSIK_IsComponentHandleValid(stiffnessWorld, tetMesh) == 1);

    PHYSIK_AddPointConnection(
        stiffnessWorld,
        node3,
        node3,
        node3,
        node3,
        1.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        2.0f,
        10.0f,
        0.0f);
    PHYSIK_Step(stiffnessWorld, 0.5f);
    const float velocityWithStiffness = GetNodeVelocity(stiffnessWorld, node3).z;

    assert(velocityWithStiffness > 0.0f);
    assert(velocityWithStiffness < velocityWithoutStiffness);

    PHYSIK_DestroyWorld(stiffnessWorld);
}

void TetMeshComponentOwnsTetsAndWorldStepUsesComponentSystem()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    const int node0 = AddNode(world, 0.0f, 0.0f, 0.0f);
    const int node1 = AddNode(world, 1.0f, 0.0f, 0.0f);
    const int node2 = AddNode(world, 0.0f, 1.0f, 0.0f);
    const int node3 = AddNode(world, 0.0f, 0.0f, 1.0f);

    const int nodes[] = {node0, node1, node2, node3};
    const int tetNodeIndices[] = {node0, node1, node2, node3};
    const PhysiK::ComponentHandle handle =
        PHYSIK_CreateTetMeshComponent(world, nodes, 4, tetNodeIndices, 1);

    assert(PHYSIK_IsComponentHandleValid(world, handle) == 1);
    assert(PHYSIK_GetTetMeshTetCount(world, handle) == 1);

    PHYSIK_Step(world, 0.01f);
    PHYSIK_DestroyWorld(world);
}

void LinearTetAssemblyProducesForcesAndSymmetricStiffness()
{
    std::vector<PhysiK::Node> nodes = CreateUnitTetNodes();
    PhysiK::Tet tet = CreateUnitTet();
    PhysiK::FEMModel::InitializeTetRestData(tet, nodes);

    std::vector<PhysiK::Tet> tets = {tet};
    PhysiK::SolverData solverData;
    PhysiK::FEMModel::AccumulateElasticForces(tets, nodes, solverData);

    for (int node = 0; node < 4; ++node)
    {
        assert(LengthSquared(SumForcesForNode(solverData, node)) < 0.000001f);
    }

    assert(solverData.GetStiffnessBlocks().size() == 16);
    for (const PhysiK::SolverData::StiffnessBlock& block : solverData.GetStiffnessBlocks())
    {
        assert(IsFinite(block.block));
    }

    for (int nodeA = 0; nodeA < 4; ++nodeA)
    {
        for (int nodeB = 0; nodeB < 4; ++nodeB)
        {
            const PhysiK::SolverData::StiffnessBlock* ab = FindBlock(solverData, nodeA, nodeB);
            const PhysiK::SolverData::StiffnessBlock* ba = FindBlock(solverData, nodeB, nodeA);
            assert(ab != nullptr);
            assert(ba != nullptr);

            for (int row = 0; row < 3; ++row)
            {
                for (int column = 0; column < 3; ++column)
                {
                    const float lhs = GetMat3Value(ab->block, row, column);
                    const float rhs = GetMat3Value(ba->block, column, row);
                    assert(std::abs(lhs - rhs) < 0.0001f);
                }
            }
        }
    }

    solverData.Clear();
    nodes[3].position = PhysiK::Vec3{0.0f, 0.0f, 1.1f};
    PhysiK::FEMModel::AccumulateElasticForces(tets, nodes, solverData);

    const PhysiK::Vec3 node3Force = SumForcesForNode(solverData, 3);
    assert(LengthSquared(node3Force) > 0.000001f);
    assert(node3Force.z < 0.0f);
}

void UnitTetShapeFunctionGradientsMatchExpectedConvention()
{
    std::vector<PhysiK::Node> nodes = CreateUnitTetNodes();
    PhysiK::Tet tet = CreateUnitTet();
    PhysiK::FEMModel::InitializeTetRestData(tet, nodes);

    assert(NearlyEqual(tet.restVolume, 1.0f / 6.0f));
    assert(NearlyEqual(tet.shapeFunctionGradients[0], PhysiK::Vec3{-1.0f, -1.0f, -1.0f}));
    assert(NearlyEqual(tet.shapeFunctionGradients[1], PhysiK::Vec3{1.0f, 0.0f, 0.0f}));
    assert(NearlyEqual(tet.shapeFunctionGradients[2], PhysiK::Vec3{0.0f, 1.0f, 0.0f}));
    assert(NearlyEqual(tet.shapeFunctionGradients[3], PhysiK::Vec3{0.0f, 0.0f, 1.0f}));
}

void DegenerateTetIsSkippedWithoutInvalidAssembly()
{
    std::vector<PhysiK::Node> nodes(4);
    nodes[0].position = PhysiK::Vec3{0.0f, 0.0f, 0.0f};
    nodes[1].position = PhysiK::Vec3{1.0f, 0.0f, 0.0f};
    nodes[2].position = PhysiK::Vec3{0.0f, 1.0f, 0.0f};
    nodes[3].position = PhysiK::Vec3{0.0f, 0.0f, 0.0f};

    PhysiK::Tet tet = CreateUnitTet();
    PhysiK::FEMModel::InitializeTetRestData(tet, nodes);

    assert(tet.restVolume == 0.0f);
    for (const PhysiK::Vec3& gradient : tet.shapeFunctionGradients)
    {
        assert(NearlyEqual(gradient, PhysiK::Vec3{}));
    }

    std::vector<PhysiK::Tet> tets = {tet};
    PhysiK::SolverData solverData;
    PhysiK::FEMModel::AccumulateElasticForces(tets, nodes, solverData);

    assert(solverData.GetNodeForces().empty());
    assert(solverData.GetStiffnessBlocks().empty());
}

void LinearTetMaterialSanitizationAvoidsInvalidForces()
{
    std::vector<PhysiK::Node> nodes = CreateUnitTetNodes();
    PhysiK::Tet noElasticity = CreateUnitTet(-10.0f, 0.25f);
    PhysiK::FEMModel::InitializeTetRestData(noElasticity, nodes);
    nodes[1].position.x += 0.1f;

    PhysiK::SolverData solverData;
    PhysiK::FEMModel::AccumulateElasticForces({noElasticity}, nodes, solverData);

    for (const PhysiK::SolverData::NodeForce& force : solverData.GetNodeForces())
    {
        assert(IsFinite(force.force));
        assert(LengthSquared(force.force) < 0.000001f);
    }

    nodes = CreateUnitTetNodes();
    PhysiK::Tet clampedPoisson = CreateUnitTet(100.0f, 0.99f);
    PhysiK::FEMModel::InitializeTetRestData(clampedPoisson, nodes);
    nodes[1].position.x += 0.1f;
    solverData.Clear();
    PhysiK::FEMModel::AccumulateElasticForces({clampedPoisson}, nodes, solverData);

    for (const PhysiK::SolverData::NodeForce& force : solverData.GetNodeForces())
    {
        assert(IsFinite(force.force));
    }

    for (const PhysiK::SolverData::StiffnessBlock& block : solverData.GetStiffnessBlocks())
    {
        assert(IsFinite(block.block));
    }
}

void LinearTetDisplacedNodeReceivesRestoringForce()
{
    std::vector<PhysiK::Node> nodes = CreateUnitTetNodes();
    PhysiK::Tet tet = CreateUnitTet();
    PhysiK::FEMModel::InitializeTetRestData(tet, nodes);
    nodes[1].position.x += 0.1f;

    PhysiK::SolverData solverData;
    PhysiK::FEMModel::AccumulateElasticForces({tet}, nodes, solverData);

    const PhysiK::Vec3 node1Force = SumForcesForNode(solverData, 1);
    assert(LengthSquared(node1Force) > 0.000001f);
    assert(node1Force.x < 0.0f);
    assert(IsFinite(node1Force));
}

void DestroyComponentInvalidatesHandle()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    int nodes[4] = {};
    CreateSingleTet(world, nodes);

    const PhysiK::ComponentHandle sphere =
        PHYSIK_CreateCollisionSphereComponent(world, 2.0f, 2.0f, 2.0f, 0.75f);
    assert(PHYSIK_IsComponentHandleValid(world, sphere) == 1);

    PHYSIK_DestroyComponent(world, sphere);
    assert(PHYSIK_IsComponentHandleValid(world, sphere) == 0);

    PHYSIK_SetCollisionComponentKinematicTarget(world, sphere, 0.25f, 0.25f, 0.25f);
    const Point before = GetTetCentroid(world, nodes);
    PHYSIK_Step(world, 0.1f);
    const Point after = GetTetCentroid(world, nodes);

    assert(DistanceSquared(after, before) < 0.000001f);

    PHYSIK_DestroyWorld(world);
}

int main()
{
    ManualPointConnectionMovesBarycentricPoint();
    GravityMovesDynamicNode();
    ImplicitEulerGravityMatchesSemiImplicitEulerForFreeNode();
    ImplicitEulerFixedNodeDoesNotMove();
    AddNodeCreatesDynamicGeometryNode();
    SetNodeFixedControlsInternalInverseMass();
    LegacyAddNodeWithInverseMassStillWorks();
    SphereContactCreatesTransientConnectionAndMovesTet();
    MultipleForceSourcesCoexist();
    ExternalLogicHookRunsOnceBeforeSubsteps();
    KinematicUpdateRunsAfterExternalLogicBeforePhysicsSubsteps();
    FEMElasticityMovesDistortedTetTowardRestShape();
    ImplicitEulerFEMTetMovesDistortedNodeTowardRestShape();
    ImplicitEulerAllDynamicTetResistsRestShapeVelocity();
    ImplicitEulerPointAnchoredTetRecoversFreeNodes();
    ImplicitEulerFEMRecoveryBeatsGravityOnlyMotion();
    FEMLumpedMassUsesDensityAndRestVolume();
    FEMLumpedMassDensityReducesAccelerationFromSameForce();
    FEMLumpedMassPreservesFixedNodes();
    FEMGravityAccelerationIsIndependentOfDensity();
    ImplicitEulerUsesStiffnessBlocks();
    TetMeshComponentOwnsTetsAndWorldStepUsesComponentSystem();
    UnitTetShapeFunctionGradientsMatchExpectedConvention();
    DegenerateTetIsSkippedWithoutInvalidAssembly();
    LinearTetMaterialSanitizationAvoidsInvalidForces();
    LinearTetDisplacedNodeReceivesRestoringForce();
    LinearTetAssemblyProducesForcesAndSymmetricStiffness();
    DestroyComponentInvalidatesHandle();
    return 0;
}
