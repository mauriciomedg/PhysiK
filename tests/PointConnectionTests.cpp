#include "PhysiK/API/PhysiKAPI.h"
#include "PhysiK/Components/SurfaceExtractionComponent.h"
#include "PhysiK/Components/SurfaceVisualComponent.h"
#include "PhysiK/Components/TetMeshComponent.h"
#include "PhysiK/Components/TetMeshPhysicsComponent.h"
#include "PhysiK/Components/TetMeshMapperComponent.h"
#include "PhysiK/Components/TopologyMeshComponent.h"
#include "PhysiK/Components/VisualMeshComponent.h"
#include "PhysiK/Geometry/TetMeshGenerator.h"
#include "PhysiK/Core/Physics/FEM/FEMModel.h"
#include "PhysiK/Core/Events/EventSystem.h"
#include "PhysiK/Core/Solvers/Linear/ConjugateGradientSolver.h"
#include "PhysiK/Core/Solvers/Linear/LinearSolver.h"
#include "PhysiK/Core/Solvers/SolverData.h"
#include "PhysiK/Core/World/World.h"
#include "PhysiK/Math/SparseBlockMatrix.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace
{
    void Require(bool condition, const char* message)
    {
        if (!condition)
        {
            std::fprintf(stderr, "PointConnectionTests failed: %s\n", message);
            std::exit(1);
        }
    }

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

    Point GetTetMeshLocalCurrentPosition(
        PhysiK::WorldHandle world,
        PhysiK::ComponentHandle tetMesh,
        int localNodeIndex)
    {
        Point point;
        const int read = PHYSIK_GetTetMeshLocalCurrentPosition(
            world,
            tetMesh,
            localNodeIndex,
            &point.x,
            &point.y,
            &point.z);
        assert(read == 1);
        return point;
    }

    bool IsNodeFixed(PhysiK::WorldHandle world, int nodeIndex)
    {
        return PHYSIK_IsNodeFixed(world, nodeIndex) != 0;
    }

    void UseGeneratedNodes(int* nodes, int nodeCount)
    {
        if (nodes == nullptr || nodeCount <= 0)
        {
            return;
        }

        const int firstGeneratedNode = nodes[nodeCount - 1] + 1;
        for (int i = 0; i < nodeCount; ++i)
        {
            nodes[i] = firstGeneratedNode + i;
        }
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

    PhysikMaterialDesc MakeMaterialDesc(
        float density,
        float youngModulus,
        float poissonRatio = 0.3f,
        float damping = 0.0f)
    {
        return PhysikMaterialDesc{density, youngModulus, poissonRatio, damping};
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

    PhysiK::Vec3 RotateZ90(const PhysiK::Vec3& value)
    {
        return PhysiK::Vec3{-value.y, value.x, value.z};
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

    PhysiK::Mat3 DiagonalBlock(float value)
    {
        return PhysiK::Mat3::FromColumns(
            PhysiK::Vec3{value, 0.0f, 0.0f},
            PhysiK::Vec3{0.0f, value, 0.0f},
            PhysiK::Vec3{0.0f, 0.0f, value});
    }

    bool HasSparseBlock(const PhysiK::SparseBlockMatrix& matrix, int rowBlock, int colBlock)
    {
        return matrix.FindBlockIndex(rowBlock, colBlock) >= 0;
    }

    std::vector<std::pair<int, int>> BuildSparsePatternFromTetConnectivity(
        const std::vector<PhysiK::Tet>& tets)
    {
        std::vector<std::pair<int, int>> blockCoordinates;
        blockCoordinates.reserve(tets.size() * 16u);

        for (const PhysiK::Tet& tet : tets)
        {
            const int nodes[4] = {tet.node0, tet.node1, tet.node2, tet.node3};
            for (int row = 0; row < 4; ++row)
            {
                for (int column = 0; column < 4; ++column)
                {
                    blockCoordinates.push_back({nodes[row], nodes[column]});
                }
            }
        }

        return blockCoordinates;
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

    float SumForceLengthSquared(const PhysiK::SolverData& solverData)
    {
        float total = 0.0f;
        for (const PhysiK::SolverData::NodeForce& force : solverData.GetNodeForces())
        {
            total += LengthSquared(force.force);
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

    std::vector<PhysiK::Vec3> RestPositionsFromNodes(
        const std::vector<PhysiK::Node>& nodes)
    {
        std::vector<PhysiK::Vec3> positions;
        positions.reserve(nodes.size());
        for (const PhysiK::Node& node : nodes)
        {
            positions.push_back(node.restPosition);
        }
        return positions;
    }

    std::vector<PhysiK::Vec3> PositionsFromNodes(
        const std::vector<PhysiK::Node>& nodes)
    {
        std::vector<PhysiK::Vec3> positions;
        positions.reserve(nodes.size());
        for (const PhysiK::Node& node : nodes)
        {
            positions.push_back(node.position);
        }
        return positions;
    }

    std::vector<PhysiK::Vec3> VelocitiesFromNodes(
        const std::vector<PhysiK::Node>& nodes)
    {
        std::vector<PhysiK::Vec3> velocities;
        velocities.reserve(nodes.size());
        for (const PhysiK::Node& node : nodes)
        {
            velocities.push_back(node.velocity);
        }
        return velocities;
    }

    void AssembleLocalFemContributions(
        const std::vector<PhysiK::TetElementContribution>& contributions,
        PhysiK::SolverData& solverData)
    {
        for (const PhysiK::TetElementContribution& contribution : contributions)
        {
            for (int node = 0; node < 4; ++node)
            {
                solverData.AddNodeForce(
                    contribution.localNodeIndices[node],
                    contribution.forces[node]);
            }

            for (int rowNode = 0; rowNode < 4; ++rowNode)
            {
                for (int columnNode = 0; columnNode < 4; ++columnNode)
                {
                    solverData.AddStiffnessBlock(
                        contribution.localNodeIndices[rowNode],
                        contribution.localNodeIndices[columnNode],
                        contribution.stiffness[rowNode][columnNode]);
                }
            }
        }
    }

    bool ComputeFemForcesIntoSolverData(
        PhysiK::FemModel femModel,
        const std::vector<PhysiK::Tet>& tets,
        const std::vector<PhysiK::Node>& nodes,
        PhysiK::SolverData& solverData)
    {
        std::vector<PhysiK::TetElementContribution> contributions;
        const bool implemented = PhysiK::FEMModel::ComputeForces(
            femModel,
            tets,
            PositionsFromNodes(nodes),
            VelocitiesFromNodes(nodes),
            contributions);
        AssembleLocalFemContributions(contributions, solverData);
        return implemented;
    }

    void ComputeLinearFemForcesIntoSolverData(
        const std::vector<PhysiK::Tet>& tets,
        const std::vector<PhysiK::Node>& nodes,
        PhysiK::SolverData& solverData)
    {
        std::vector<PhysiK::TetElementContribution> contributions;
        PhysiK::FEMModel::ComputeElasticForces(
            tets,
            PositionsFromNodes(nodes),
            VelocitiesFromNodes(nodes),
            contributions);
        AssembleLocalFemContributions(contributions, solverData);
    }

    std::vector<PhysiK::Node> CreateUnitTetNodes()
    {
        std::vector<PhysiK::Node> nodes(4);
        nodes[0].position = PhysiK::Vec3{0.0f, 0.0f, 0.0f};
        nodes[1].position = PhysiK::Vec3{1.0f, 0.0f, 0.0f};
        nodes[2].position = PhysiK::Vec3{0.0f, 1.0f, 0.0f};
        nodes[3].position = PhysiK::Vec3{0.0f, 0.0f, 1.0f};
        for (PhysiK::Node& node : nodes)
        {
            node.restPosition = node.position;
        }
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

    std::vector<PhysiK::Node> CreateTwoTetNodes()
    {
        std::vector<PhysiK::Node> nodes(5);
        nodes[0].position = PhysiK::Vec3{0.0f, 0.0f, 0.0f};
        nodes[1].position = PhysiK::Vec3{1.0f, 0.0f, 0.0f};
        nodes[2].position = PhysiK::Vec3{0.0f, 1.0f, 0.0f};
        nodes[3].position = PhysiK::Vec3{0.0f, 0.0f, 1.0f};
        nodes[4].position = PhysiK::Vec3{0.0f, 0.0f, -1.0f};
        for (PhysiK::Node& node : nodes)
        {
            node.restPosition = node.position;
        }
        return nodes;
    }

    PhysiK::Tet CreateLowerUnitTet(float youngModulus = 100.0f)
    {
        PhysiK::Tet tet;
        tet.node0 = 0;
        tet.node1 = 1;
        tet.node2 = 2;
        tet.node3 = 4;
        tet.youngModulus = youngModulus;
        tet.poissonRatio = 0.25f;
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
        PhysikMaterialDesc material = MakeMaterialDesc(24.0f, 25.0f, 0.3f, 0.25f);
        const PhysiK::ComponentHandle tetMesh =
            PHYSIK_CreateTetMeshPhysicsComponent(
                world,
                outNodes,
                4,
                tetNodeIndices,
                1,
                &material,
                0);
        assert(PHYSIK_IsComponentHandleValid(world, tetMesh) == 1);

        UseGeneratedNodes(outNodes, 4);
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
        PhysikMaterialDesc material = MakeMaterialDesc(density, youngModulus, 0.3f, damping);
        const PhysiK::ComponentHandle tetMesh =
            PHYSIK_CreateTetMeshPhysicsComponent(
                world,
                outNodes,
                4,
                tetNodeIndices,
                1,
                &material,
                0);
        assert(PHYSIK_IsComponentHandleValid(world, tetMesh) == 1);

        UseGeneratedNodes(outNodes, 4);
    }

    PhysiK::ComponentHandle CreateSingleTetMesh(
        PhysiK::WorldHandle world,
        int (&outNodes)[4])
    {
        outNodes[0] = AddNode(world, 0.0f, 0.0f, 0.0f);
        outNodes[1] = AddNode(world, 1.0f, 0.0f, 0.0f);
        outNodes[2] = AddNode(world, 0.0f, 1.0f, 0.0f);
        outNodes[3] = AddNode(world, 0.0f, 0.0f, 1.0f);

        const int tetNodeIndices[] = {outNodes[0], outNodes[1], outNodes[2], outNodes[3]};
        PhysikMaterialDesc material = MakeMaterialDesc(1.0f, 25.0f, 0.3f, 0.0f);
        const PhysiK::ComponentHandle tetMesh =
            PHYSIK_CreateTetMeshPhysicsComponent(
                world,
                outNodes,
                4,
                tetNodeIndices,
                1,
                &material,
                0);
        assert(PHYSIK_IsComponentHandleValid(world, tetMesh) == 1);
        UseGeneratedNodes(outNodes, 4);
        return tetMesh;
    }

    PhysiK::ComponentHandle CreateTwoTetMesh(
        PhysiK::WorldHandle world,
        int (&outNodes)[5],
        float youngModulus = 0.0f,
        float density = 1.0f)
    {
        outNodes[0] = AddNode(world, 0.0f, 0.0f, 0.0f);
        outNodes[1] = AddNode(world, 1.0f, 0.0f, 0.0f);
        outNodes[2] = AddNode(world, 0.0f, 1.0f, 0.0f);
        outNodes[3] = AddNode(world, 0.0f, 0.0f, 1.0f);
        outNodes[4] = AddNode(world, 0.0f, 0.0f, -1.0f);

        const int tetNodeIndices[] = {
            outNodes[0], outNodes[1], outNodes[2], outNodes[3],
            outNodes[0], outNodes[1], outNodes[2], outNodes[4]};
        PhysikMaterialDesc material = MakeMaterialDesc(density, youngModulus, 0.3f, 0.0f);
        const PhysiK::ComponentHandle tetMesh =
            PHYSIK_CreateTetMeshPhysicsComponent(
                world,
                outNodes,
                5,
                tetNodeIndices,
                2,
                &material,
                0);
        assert(PHYSIK_IsComponentHandleValid(world, tetMesh) == 1);
        UseGeneratedNodes(outNodes, 5);
        return tetMesh;
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

    class RecordingEventComponent final : public PhysiK::Component
    {
    public:
        int eventCount = 0;
        PhysiK::PhysicsEvent lastEvent{
            PhysiK::PhysicsEventType::TetMeshTopologyChanged,
            nullptr,
            nullptr};

        void OnPhysicsEvent(const PhysiK::PhysicsEvent& event) override
        {
            ++eventCount;
            lastEvent = event;
        }
    };

    bool TriangleWindsAwayFromOppositePoint(
        const PhysiK::Vec3* positions,
        int node0,
        int node1,
        int node2,
        int oppositeNode)
    {
        const PhysiK::Vec3& p0 = positions[node0];
        const PhysiK::Vec3& p1 = positions[node1];
        const PhysiK::Vec3& p2 = positions[node2];
        const PhysiK::Vec3& po = positions[oppositeNode];
        const PhysiK::Vec3 normal = PhysiK::Cross(p1 - p0, p2 - p0);
        return PhysiK::Dot(normal, po - p0) <= 0.0f;
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

    PHYSIK_DestroyWorld(world);
}

void SetNodeFixedControlsDynamicState()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    const int node = PHYSIK_AddNode(world, 0.0f, 0.0f, 0.0f);

    PHYSIK_SetNodeFixed(world, node, 1);
    assert(IsNodeFixed(world, node));

    PHYSIK_SetNodeFixed(world, node, 0);
    assert(!IsNodeFixed(world, node));

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

void CollisionSphereConnectionSettingsHaveDefaults()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    const PhysiK::ComponentHandle sphere =
        PHYSIK_CreateCollisionSphereComponent(world, 0.0f, 0.0f, 0.0f, 0.5f);

    float stiffness = -1.0f;
    float damping = -1.0f;
    PHYSIK_GetCollisionSphereConnectionSettings(world, sphere, &stiffness, &damping);

    assert(NearlyEqual(stiffness, 1000.0f));
    assert(NearlyEqual(damping, 10.0f));

    PHYSIK_DestroyWorld(world);
}

void CollisionSphereConnectionSettingsCanBeUpdatedAndClamped()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    const PhysiK::ComponentHandle sphere =
        PHYSIK_CreateCollisionSphereComponent(world, 0.0f, 0.0f, 0.0f, 0.5f);

    PHYSIK_SetCollisionSphereConnectionSettings(world, sphere, 2500.0f, 12.5f);

    float stiffness = 0.0f;
    float damping = 0.0f;
    PHYSIK_GetCollisionSphereConnectionSettings(world, sphere, &stiffness, &damping);
    assert(NearlyEqual(stiffness, 2500.0f));
    assert(NearlyEqual(damping, 12.5f));

    PHYSIK_SetCollisionSphereConnectionSettings(world, sphere, -10.0f, -2.0f);
    PHYSIK_GetCollisionSphereConnectionSettings(world, sphere, &stiffness, &damping);
    assert(NearlyEqual(stiffness, 0.0f));
    assert(NearlyEqual(damping, 0.0f));

    PHYSIK_DestroyWorld(world);
}

void CollisionSphereConnectionSettingsAffectGeneratedConnections()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    int nodes[4] = {};
    CreateSingleTet(world, nodes);
    const PhysiK::ComponentHandle sphere =
        PHYSIK_CreateCollisionSphereComponent(world, 0.25f, 0.25f, 0.25f, 0.75f);
    PHYSIK_SetCollisionSphereConnectionSettings(world, sphere, 0.0f, 0.0f);

    const Point before = GetTetCentroid(world, nodes);
    PHYSIK_Step(world, 0.1f);
    const Point after = GetTetCentroid(world, nodes);

    assert(DistanceSquared(before, after) < 0.000001f);

    PHYSIK_DestroyWorld(world);
}

void CollisionSphereConnectionSettingsCAPIHandlesInvalidInputs()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    int nodes[4] = {};
    const PhysiK::ComponentHandle tetMesh = CreateSingleTetMesh(world, nodes);
    const PhysiK::ComponentHandle sphere =
        PHYSIK_CreateCollisionSphereComponent(world, 0.0f, 0.0f, 0.0f, 0.5f);

    float stiffness = -1.0f;
    float damping = -1.0f;

    PHYSIK_SetCollisionSphereConnectionSettings(nullptr, sphere, 10.0f, 1.0f);
    PHYSIK_SetCollisionSphereConnectionSettings(world, PhysiK::ComponentHandle{}, 10.0f, 1.0f);
    PHYSIK_SetCollisionSphereConnectionSettings(world, tetMesh, 10.0f, 1.0f);

    PHYSIK_GetCollisionSphereConnectionSettings(nullptr, sphere, &stiffness, &damping);
    assert(NearlyEqual(stiffness, -1.0f));
    assert(NearlyEqual(damping, -1.0f));

    PHYSIK_GetCollisionSphereConnectionSettings(world, tetMesh, &stiffness, &damping);
    assert(NearlyEqual(stiffness, -1.0f));
    assert(NearlyEqual(damping, -1.0f));

    PHYSIK_GetCollisionSphereConnectionSettings(world, sphere, &stiffness, nullptr);
    assert(NearlyEqual(stiffness, 1000.0f));

    PHYSIK_GetCollisionSphereConnectionSettings(world, sphere, nullptr, &damping);
    assert(NearlyEqual(damping, 10.0f));

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

    int nodes[] = {node0, node1, node2, node3};
    const int tetNodeIndices[] = {node0, node1, node2, node3};
    PhysikMaterialDesc material = MakeMaterialDesc(24.0f, 25.0f, 0.3f, 0.25f);
    const PhysiK::ComponentHandle tetMesh =
        PHYSIK_CreateTetMeshPhysicsComponent(
            world,
            nodes,
            4,
            tetNodeIndices,
            1,
            &material,
            0);
    assert(PHYSIK_IsComponentHandleValid(world, tetMesh) == 1);
    UseGeneratedNodes(nodes, 4);
    PHYSIK_SetNodeFixed(world, nodes[0], 1);
    PHYSIK_SetNodeFixed(world, nodes[1], 1);
    PHYSIK_SetNodeFixed(world, nodes[2], 1);

    const Point restPosition = GetNodePosition(world, nodes[3]);
    PHYSIK_SetNodePosition(world, nodes[3], 0.0f, 0.0f, 1.25f);
    const Point distortedPosition = GetNodePosition(world, nodes[3]);

    PHYSIK_SetSubstepCount(world, 4);
    PHYSIK_Step(world, 0.1f);

    const Point after = GetNodePosition(world, nodes[3]);

    assert(DistanceSquared(after, restPosition) < DistanceSquared(distortedPosition, restPosition));
    assert(after.z < distortedPosition.z);

    PHYSIK_DestroyWorld(world);
}

void FEMUsesWorldNodeMappingWhenLocalTetIndicesDiffer()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);
    PHYSIK_SetGravity(world, 0.0f, 0.0f, 0.0f);

    AddNode(world, 100.0f, 100.0f, 100.0f);
    const int node0 = AddFixedNode(world, 0.0f, 0.0f, 0.0f);
    const int node1 = AddFixedNode(world, 1.0f, 0.0f, 0.0f);
    const int node2 = AddFixedNode(world, 0.0f, 1.0f, 0.0f);
    const int node3 = AddNode(world, 0.0f, 0.0f, 1.0f);

    int nodes[] = {node0, node1, node2, node3};
    const int tetNodeIndices[] = {node0, node1, node2, node3};
    PhysikMaterialDesc material = MakeMaterialDesc(24.0f, 25.0f, 0.3f, 0.25f);
    const PhysiK::ComponentHandle tetMesh =
        PHYSIK_CreateTetMeshPhysicsComponent(
            world,
            nodes,
            4,
            tetNodeIndices,
            1,
            &material,
            0);
    assert(PHYSIK_IsComponentHandleValid(world, tetMesh) == 1);
    assert(node0 == 1);
    assert(node3 == 4);
    UseGeneratedNodes(nodes, 4);
    PHYSIK_SetNodeFixed(world, nodes[0], 1);
    PHYSIK_SetNodeFixed(world, nodes[1], 1);
    PHYSIK_SetNodeFixed(world, nodes[2], 1);

    PHYSIK_SetNodePosition(world, nodes[3], 0.0f, 0.0f, 1.25f);
    PHYSIK_Step(world, 0.01f);

    const Point after = GetNodePosition(world, nodes[3]);
    assert(after.z < 1.25f);

    const Point dummy = GetNodePosition(world, 0);
    assert(NearlyEqual(dummy.x, 100.0f));
    assert(NearlyEqual(dummy.y, 100.0f));
    assert(NearlyEqual(dummy.z, 100.0f));

    PHYSIK_DestroyWorld(world);
}

void TetMeshPhysicsUsesCreationRestDataForRuntimeForces()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);
    PHYSIK_SetGravity(world, 0.0f, 0.0f, 0.0f);

    const int node0 = AddFixedNode(world, 0.0f, 0.0f, 0.0f);
    const int node1 = AddFixedNode(world, 1.0f, 0.0f, 0.0f);
    const int node2 = AddFixedNode(world, 0.0f, 1.0f, 0.0f);
    const int node3 = AddNode(world, 0.0f, 0.0f, 1.0f);

    int nodes[] = {node0, node1, node2, node3};
    const int tetNodeIndices[] = {node0, node1, node2, node3};
    PhysikMaterialDesc material = MakeMaterialDesc(24.0f, 25.0f, 0.3f, 0.0f);
    const PhysiK::ComponentHandle tetMesh =
        PHYSIK_CreateTetMeshPhysicsComponent(
            world,
            nodes,
            4,
            tetNodeIndices,
            1,
            &material,
            0);
    assert(PHYSIK_IsComponentHandleValid(world, tetMesh) == 1);
    UseGeneratedNodes(nodes, 4);
    PHYSIK_SetNodeFixed(world, nodes[0], 1);
    PHYSIK_SetNodeFixed(world, nodes[1], 1);
    PHYSIK_SetNodeFixed(world, nodes[2], 1);

    PHYSIK_SetNodePosition(world, nodes[3], 0.0f, 0.0f, 1.25f);
    PHYSIK_Step(world, 0.01f);
    const float firstStepZ = GetNodePosition(world, nodes[3]).z;
    assert(firstStepZ < 1.25f);

    PHYSIK_SetNodePosition(world, nodes[3], 0.0f, 0.0f, 1.50f);
    PHYSIK_SetNodeVelocity(world, nodes[3], 0.0f, 0.0f, 0.0f);
    PHYSIK_Step(world, 0.01f);
    const float secondStepZ = GetNodePosition(world, nodes[3]).z;
    assert(secondStepZ < 1.50f);

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

    int nodes[] = {node0, node1, node2, node3};
    const int tetNodeIndices[] = {node0, node1, node2, node3};
    PhysikMaterialDesc material = MakeMaterialDesc(24.0f, 25.0f, 0.3f, 0.25f);
    const PhysiK::ComponentHandle tetMesh =
        PHYSIK_CreateTetMeshPhysicsComponent(
            world,
            nodes,
            4,
            tetNodeIndices,
            1,
            &material,
            0);
    assert(PHYSIK_IsComponentHandleValid(world, tetMesh) == 1);
    UseGeneratedNodes(nodes, 4);
    PHYSIK_SetNodeFixed(world, nodes[0], 1);
    PHYSIK_SetNodeFixed(world, nodes[1], 1);
    PHYSIK_SetNodeFixed(world, nodes[2], 1);

    const Point restPosition = GetNodePosition(world, nodes[3]);
    PHYSIK_SetNodePosition(world, nodes[3], 0.0f, 0.0f, 1.25f);
    const Point distortedPosition = GetNodePosition(world, nodes[3]);

    PHYSIK_SetSolverMode(world, 1);
    PHYSIK_Step(world, 0.25f);

    const Point after = GetNodePosition(world, nodes[3]);

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

    const Point velocity = GetNodeVelocity(world, nodes[3]);
    assert(NearlyEqual(velocity.z, 2.4f, 0.0001f));

    PHYSIK_DestroyWorld(world);
}

void TetMeshMaterialDescCrossesNativeBoundary()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    const int node0 = AddNode(world, 0.0f, 0.0f, 0.0f);
    const int node1 = AddNode(world, 1.0f, 0.0f, 0.0f);
    const int node2 = AddNode(world, 0.0f, 1.0f, 0.0f);
    const int node3 = AddNode(world, 0.0f, 0.0f, 1.0f);
    int nodes[] = {node0, node1, node2, node3};
    const int tetNodeIndices[] = {node0, node1, node2, node3};
    PhysikMaterialDesc material = MakeMaterialDesc(1.0f, 0.0f);

    const PhysiK::ComponentHandle invalidTetMesh =
        PHYSIK_CreateTetMeshPhysicsComponent(
            world,
            nodes,
            4,
            tetNodeIndices,
            1,
            nullptr,
            0);
    assert(PHYSIK_IsComponentHandleValid(world, invalidTetMesh) == 0);

    const PhysiK::ComponentHandle tetMesh =
        PHYSIK_CreateTetMeshPhysicsComponent(
            world,
            nodes,
            4,
            tetNodeIndices,
            1,
            &material,
            0);
    assert(PHYSIK_IsComponentHandleValid(world, tetMesh) == 1);
    UseGeneratedNodes(nodes, 4);

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

    const Point velocity = GetNodeVelocity(world, nodes[3]);
    assert(NearlyEqual(velocity.z, 2.4f, 0.0001f));

    PHYSIK_DestroyWorld(world);
}

void TetMeshMaterialCanBeUpdatedThroughNativeDescriptor()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    const int node0 = AddNode(world, 0.0f, 0.0f, 0.0f);
    const int node1 = AddNode(world, 1.0f, 0.0f, 0.0f);
    const int node2 = AddNode(world, 0.0f, 1.0f, 0.0f);
    const int node3 = AddNode(world, 0.0f, 0.0f, 1.0f);
    int nodes[] = {node0, node1, node2, node3};
    const int tetNodeIndices[] = {node0, node1, node2, node3};
    PhysikMaterialDesc material = MakeMaterialDesc(1.0f, 0.0f);

    const PhysiK::ComponentHandle tetMesh =
        PHYSIK_CreateTetMeshPhysicsComponent(
            world,
            nodes,
            4,
            tetNodeIndices,
            1,
            &material,
            0);
    assert(PHYSIK_IsComponentHandleValid(world, tetMesh) == 1);
    UseGeneratedNodes(nodes, 4);

    PhysikMaterialDesc heavierMaterial = MakeMaterialDesc(4.0f, 0.0f);
    PHYSIK_SetTetMeshMaterial(world, tetMesh, nullptr);
    PHYSIK_SetTetMeshMaterial(world, tetMesh, &heavierMaterial);

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

    const Point velocity = GetNodeVelocity(world, nodes[3]);
    assert(NearlyEqual(velocity.z, 0.6f, 0.0001f));

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
    int nodes[] = {node0, node1, node2, node3};
    const int tetNodeIndices[] = {node0, node1, node2, node3};
    PhysikMaterialDesc material = MakeMaterialDesc(1.0f, 0.0f, 0.3f, 0.0f);

    const PhysiK::ComponentHandle tetMesh =
        PHYSIK_CreateTetMeshPhysicsComponent(
            world,
            nodes,
            4,
            tetNodeIndices,
            1,
            &material,
            0);
    assert(PHYSIK_IsComponentHandleValid(world, tetMesh) == 1);
    UseGeneratedNodes(nodes, 4);

    assert(!IsNodeFixed(world, nodes[0]));
    assert(!IsNodeFixed(world, nodes[1]));
    assert(!IsNodeFixed(world, nodes[2]));
    assert(!IsNodeFixed(world, nodes[3]));

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
    int componentNodes[] = {node0, node1, node2, node3};
    const int tetNodeIndices[] = {node0, node1, node2, node3};
    PhysikMaterialDesc material = MakeMaterialDesc(24.0f, 100.0f, 0.3f, 0.0f);
    const PhysiK::ComponentHandle tetMesh =
        PHYSIK_CreateTetMeshPhysicsComponent(
            stiffnessWorld,
            componentNodes,
            4,
            tetNodeIndices,
            1,
            &material,
            0);
    assert(PHYSIK_IsComponentHandleValid(stiffnessWorld, tetMesh) == 1);
    UseGeneratedNodes(componentNodes, 4);
    PHYSIK_SetNodeFixed(stiffnessWorld, componentNodes[0], 1);
    PHYSIK_SetNodeFixed(stiffnessWorld, componentNodes[1], 1);
    PHYSIK_SetNodeFixed(stiffnessWorld, componentNodes[2], 1);

    PHYSIK_AddPointConnection(
        stiffnessWorld,
        componentNodes[3],
        componentNodes[3],
        componentNodes[3],
        componentNodes[3],
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
    const float velocityWithStiffness =
        GetNodeVelocity(stiffnessWorld, componentNodes[3]).z;

    assert(velocityWithStiffness > 0.0f);
    assert(velocityWithStiffness < velocityWithoutStiffness);

    PHYSIK_DestroyWorld(stiffnessWorld);
}

void ImplicitPointConnectionStableWithHighStiffness()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    PHYSIK_SetSolverMode(world, 1);
    PHYSIK_SetGravity(world, 0.0f, 0.0f, 0.0f);
    const int node = AddNode(world, 0.0f, 0.0f, 0.0f);

    PHYSIK_AddPointConnection(
        world,
        node,
        node,
        node,
        node,
        1.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        1.0f,
        10000.0f,
        0.0f);
    assert(PHYSIK_GetPointConnectionCount(world) == 1);

    PHYSIK_Step(world, 0.1f);

    const Point position = GetNodePosition(world, node);
    const Point velocity = GetNodeVelocity(world, node);
    assert(IsFinite(position.z));
    assert(IsFinite(velocity.z));
    assert(position.z > 0.9f);
    assert(position.z < 1.1f);
    assert(PHYSIK_GetPointConnectionCount(world) == 0);

    PHYSIK_DestroyWorld(world);
}

void PointConnectionBarycentricAssemblyDistributesForcesAndStiffness()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    int nodes[4] = {};
    CreateSingleTetWithMaterial(world, nodes, 0.0f, 0.0f, 24.0f);
    PHYSIK_SetSolverMode(world, 1);
    PHYSIK_SetGravity(world, 0.0f, 0.0f, 0.0f);

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

    PHYSIK_Step(world, 0.1f);

    const float velocity0 = GetNodeVelocity(world, nodes[0]).z;
    const float velocity1 = GetNodeVelocity(world, nodes[1]).z;
    const float velocity2 = GetNodeVelocity(world, nodes[2]).z;
    const float velocity3 = GetNodeVelocity(world, nodes[3]).z;

    assert(velocity0 > 0.0f);
    assert(NearlyEqual(velocity0, velocity1, 0.0001f));
    assert(NearlyEqual(velocity0, velocity2, 0.0001f));
    assert(NearlyEqual(velocity0, velocity3, 0.0001f));

    PHYSIK_DestroyWorld(world);
}

void ImplicitAnchoredTetPointConnectionsRemainStableUnderGravity()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    int nodes[4] = {};
    CreateSingleTetWithMaterial(world, nodes, 250.0f, 0.0f, 1.0f);
    PHYSIK_SetSolverMode(world, 1);
    PHYSIK_SetGravity(world, 0.0f, -10.0f, 0.0f);

    const Point anchors[3] = {
        GetNodePosition(world, nodes[0]),
        GetNodePosition(world, nodes[1]),
        GetNodePosition(world, nodes[2])};

    for (int i = 0; i < 3; ++i)
    {
        PHYSIK_AddPointConnection(
            world,
            nodes[i],
            nodes[i],
            nodes[i],
            nodes[i],
            1.0f,
            0.0f,
            0.0f,
            0.0f,
            anchors[i].x,
            anchors[i].y,
            anchors[i].z,
            10000.0f,
            0.0f);
    }

    PHYSIK_Step(world, 0.1f);

    for (int i = 0; i < 4; ++i)
    {
        const Point position = GetNodePosition(world, nodes[i]);
        const Point velocity = GetNodeVelocity(world, nodes[i]);
        assert(IsFinite(position.x));
        assert(IsFinite(position.y));
        assert(IsFinite(position.z));
        assert(IsFinite(velocity.x));
        assert(IsFinite(velocity.y));
        assert(IsFinite(velocity.z));
    }

    assert(DistanceSquared(GetNodePosition(world, nodes[0]), anchors[0]) < 0.02f);
    assert(DistanceSquared(GetNodePosition(world, nodes[1]), anchors[1]) < 0.02f);
    assert(DistanceSquared(GetNodePosition(world, nodes[2]), anchors[2]) < 0.02f);
    assert(PHYSIK_GetPointConnectionCount(world) == 0);

    PHYSIK_DestroyWorld(world);
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
    PhysikMaterialDesc material = MakeMaterialDesc(1.0f, 25.0f, 0.3f, 0.25f);
    const PhysiK::ComponentHandle handle =
        PHYSIK_CreateTetMeshPhysicsComponent(world, nodes, 4, tetNodeIndices, 1, &material, 0);

    assert(PHYSIK_IsComponentHandleValid(world, handle) == 1);
    assert(PHYSIK_GetTetMeshTetCount(world, handle) == 1);

    PHYSIK_Step(world, 0.01f);
    PHYSIK_DestroyWorld(world);
}

void CollisionSphereOverlapQueryOverlapsOneTetNode()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    int nodes[4] = {};
    const PhysiK::ComponentHandle tetMesh = CreateSingleTetMesh(world, nodes);
    const PhysiK::ComponentHandle sphere =
        PHYSIK_CreateCollisionSphereComponent(world, 0.0f, 0.0f, 0.0f, 0.05f);

    assert(PHYSIK_GetCollisionSphereOverlapCount(world, sphere) == 1);
    PhysikCollisionSphereOverlap overlaps[1] = {};
    assert(PHYSIK_GetCollisionSphereOverlaps(world, sphere, overlaps, 1) == 1);
    assert(overlaps[0].geometryType == PHYSIK_OverlapGeometry_Tetrahedron);
    assert(overlaps[0].component.index == tetMesh.index);
    assert(overlaps[0].component.generation == tetMesh.generation);
    assert(overlaps[0].primitiveIndex == 0);
    assert(overlaps[0].node0 == nodes[0]);
    assert(overlaps[0].node1 == nodes[1]);
    assert(overlaps[0].node2 == nodes[2]);
    assert(overlaps[0].node3 == nodes[3]);
    assert(overlaps[0].overlappedNodeMask == 1);
    assert(overlaps[0].overlappedNodeCount == 1);
    assert(NearlyEqual(overlaps[0].sphereCenterX, 0.0f));
    assert(NearlyEqual(overlaps[0].sphereRadius, 0.05f));
    assert(NearlyEqual(overlaps[0].minDistance, 0.0f));

    PHYSIK_DestroyWorld(world);
}

void CollisionSphereOverlapUsesWorldNodeMapping()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    AddNode(world, 100.0f, 100.0f, 100.0f);
    int nodes[4] = {};
    const PhysiK::ComponentHandle tetMesh = CreateSingleTetMesh(world, nodes);
    const PhysiK::ComponentHandle sphere =
        PHYSIK_CreateCollisionSphereComponent(world, 0.0f, 0.0f, 0.0f, 0.05f);

    assert(nodes[0] != 0);
    assert(PHYSIK_GetCollisionSphereOverlapCount(world, sphere) == 1);
    PhysikCollisionSphereOverlap overlaps[1] = {};
    assert(PHYSIK_GetCollisionSphereOverlaps(world, sphere, overlaps, 1) == 1);
    assert(overlaps[0].component.index == tetMesh.index);
    assert(overlaps[0].node0 == nodes[0]);
    assert(overlaps[0].node1 == nodes[1]);
    assert(overlaps[0].node2 == nodes[2]);
    assert(overlaps[0].node3 == nodes[3]);
    assert(overlaps[0].node0 != 0);

    PHYSIK_DestroyWorld(world);
}

void CollisionSphereOverlapQueryReturnsZeroWhenNoNodesOverlap()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    int nodes[4] = {};
    CreateSingleTetMesh(world, nodes);
    const PhysiK::ComponentHandle sphere =
        PHYSIK_CreateCollisionSphereComponent(world, 10.0f, 10.0f, 10.0f, 0.5f);

    assert(PHYSIK_GetCollisionSphereOverlapCount(world, sphere) == 0);
    PhysikCollisionSphereOverlap overlaps[1] = {};
    assert(PHYSIK_GetCollisionSphereOverlaps(world, sphere, overlaps, 1) == 0);

    PHYSIK_DestroyWorld(world);
}

void CollisionSphereOverlapQueryIgnoresInactiveTets()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    int nodes[5] = {};
    const PhysiK::ComponentHandle tetMesh = CreateTwoTetMesh(world, nodes);
    const PhysiK::ComponentHandle sphere =
        PHYSIK_CreateCollisionSphereComponent(world, 0.0f, 0.0f, -1.0f, 0.05f);

    assert(PHYSIK_GetCollisionSphereOverlapCount(world, sphere) == 1);
    PHYSIK_DeactivateTet(world, tetMesh, 1);
    assert(PHYSIK_GetCollisionSphereOverlapCount(world, sphere) == 0);

    PHYSIK_DestroyWorld(world);
}

void CollisionSphereOverlapQueryIgnoresDestroyedTetMesh()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    int nodes[4] = {};
    const PhysiK::ComponentHandle tetMesh = CreateSingleTetMesh(world, nodes);
    const PhysiK::ComponentHandle sphere =
        PHYSIK_CreateCollisionSphereComponent(world, 0.0f, 0.0f, 0.0f, 0.05f);

    assert(PHYSIK_GetCollisionSphereOverlapCount(world, sphere) == 1);
    PHYSIK_DestroyComponent(world, tetMesh);
    assert(PHYSIK_GetCollisionSphereOverlapCount(world, sphere) == 0);

    PHYSIK_DestroyWorld(world);
}

void CollisionSphereOverlapQueryIgnoresZeroRadiusSphere()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    int nodes[4] = {};
    CreateSingleTetMesh(world, nodes);
    const PhysiK::ComponentHandle sphere =
        PHYSIK_CreateCollisionSphereComponent(world, 0.0f, 0.0f, 0.0f, 0.0f);

    assert(PHYSIK_GetCollisionSphereOverlapCount(world, sphere) == 0);

    PHYSIK_DestroyWorld(world);
}

void CollisionSphereOverlapQueryReportsMaskAndCount()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    int nodes[4] = {};
    CreateSingleTetMesh(world, nodes);
    const PhysiK::ComponentHandle sphere =
        PHYSIK_CreateCollisionSphereComponent(world, 0.5f, 0.0f, 0.0f, 0.51f);

    PhysikCollisionSphereOverlap overlaps[1] = {};
    assert(PHYSIK_GetCollisionSphereOverlaps(world, sphere, overlaps, 1) == 1);
    assert(overlaps[0].overlappedNodeMask == 3);
    assert(overlaps[0].overlappedNodeCount == 2);
    assert(NearlyEqual(overlaps[0].minDistance, 0.5f));

    PHYSIK_DestroyWorld(world);
}

void CollisionSphereOverlapQueryCountAndFillAreConsistent()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    int nodes[5] = {};
    CreateTwoTetMesh(world, nodes);
    const PhysiK::ComponentHandle sphere =
        PHYSIK_CreateCollisionSphereComponent(world, 0.0f, 0.0f, 0.0f, 0.05f);

    assert(PHYSIK_GetCollisionSphereOverlapCount(world, sphere) == 2);
    PhysikCollisionSphereOverlap oneOverlap[1] = {};
    assert(PHYSIK_GetCollisionSphereOverlaps(world, sphere, oneOverlap, 1) == 1);
    PhysikCollisionSphereOverlap overlaps[2] = {};
    assert(PHYSIK_GetCollisionSphereOverlaps(world, sphere, overlaps, 2) == 2);
    assert(overlaps[0].primitiveIndex != overlaps[1].primitiveIndex);

    PHYSIK_DestroyWorld(world);
}

void CollisionSphereOverlapQueryInvalidInputsReturnZero()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    int nodes[4] = {};
    const PhysiK::ComponentHandle tetMesh = CreateSingleTetMesh(world, nodes);
    const PhysiK::ComponentHandle sphere =
        PHYSIK_CreateCollisionSphereComponent(world, 0.0f, 0.0f, 0.0f, 0.05f);
    PhysikCollisionSphereOverlap overlaps[1] = {};

    assert(PHYSIK_GetCollisionSphereOverlapCount(nullptr, sphere) == 0);
    assert(PHYSIK_GetCollisionSphereOverlapCount(world, PhysiK::ComponentHandle{}) == 0);
    assert(PHYSIK_GetCollisionSphereOverlapCount(world, tetMesh) == 0);
    assert(PHYSIK_GetCollisionSphereOverlaps(nullptr, sphere, overlaps, 1) == 0);
    assert(PHYSIK_GetCollisionSphereOverlaps(world, sphere, nullptr, 1) == 0);
    assert(PHYSIK_GetCollisionSphereOverlaps(world, sphere, overlaps, 0) == 0);
    assert(PHYSIK_GetCollisionSphereOverlaps(world, tetMesh, overlaps, 1) == 0);

    PHYSIK_DestroyWorld(world);
}

void SparseBlockMatrixStoresAndMultipliesBlocks()
{
    PhysiK::SparseBlockMatrix matrix;
    matrix.BuildPattern(2, {{0, 0}, {0, 1}, {1, 1}});
    assert(matrix.blockCount == 2);
    assert(matrix.values.size() == 3);

    assert(matrix.AddBlock(0, 0, DiagonalBlock(2.0f)));
    assert(matrix.AddBlock(0, 1, DiagonalBlock(3.0f)));
    assert(matrix.AddBlock(1, 1, DiagonalBlock(4.0f)));
    assert(!matrix.AddBlock(1, 0, DiagonalBlock(5.0f)));

    std::vector<float> input = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f};
    std::vector<float> output;
    matrix.Multiply(input, output);

    assert(output.size() == 6);
    assert(NearlyEqual(output[0], 14.0f));
    assert(NearlyEqual(output[1], 19.0f));
    assert(NearlyEqual(output[2], 24.0f));
    assert(NearlyEqual(output[3], 16.0f));
    assert(NearlyEqual(output[4], 20.0f));
    assert(NearlyEqual(output[5], 24.0f));
}

void SparseBlockMatrixAddBlockAccumulatesContributions()
{
    PhysiK::SparseBlockMatrix matrix;
    matrix.BuildPattern(1, {{0, 0}});
    assert(matrix.AddBlock(0, 0, DiagonalBlock(2.0f)));
    assert(matrix.AddBlock(0, 0, DiagonalBlock(3.0f)));

    const int blockIndex = matrix.FindBlockIndex(0, 0);
    assert(blockIndex >= 0);
    assert(NearlyEqual(GetMat3Value(matrix.values[static_cast<std::size_t>(blockIndex)], 0, 0), 5.0f));
    assert(NearlyEqual(GetMat3Value(matrix.values[static_cast<std::size_t>(blockIndex)], 1, 1), 5.0f));
    assert(NearlyEqual(GetMat3Value(matrix.values[static_cast<std::size_t>(blockIndex)], 2, 2), 5.0f));

    matrix.ClearValues();
    assert(NearlyEqual(GetMat3Value(matrix.values[static_cast<std::size_t>(blockIndex)], 0, 0), 0.0f));

    assert(matrix.AddBlock(0, 0, DiagonalBlock(7.0f)));
    assert(NearlyEqual(GetMat3Value(matrix.values[static_cast<std::size_t>(blockIndex)], 2, 2), 7.0f));
}

void SparseBlockMatrixSingleTetPatternContainsAllCouplings()
{
    PhysiK::SparseBlockMatrix matrix;
    PhysiK::Tet tet = CreateUnitTet();
    matrix.BuildPattern(4, BuildSparsePatternFromTetConnectivity({tet}));

    assert(matrix.blockCount == 4);
    assert(matrix.values.size() == 16);
    for (int row = 0; row < 4; ++row)
    {
        assert(matrix.rowStart[static_cast<std::size_t>(row + 1)] -
            matrix.rowStart[static_cast<std::size_t>(row)] == 4);
        for (int column = 0; column < 4; ++column)
        {
            assert(HasSparseBlock(matrix, row, column));
        }
    }
}

void SparseBlockMatrixAdjacentTetsReuseSharedBlocks()
{
    PhysiK::SparseBlockMatrix matrix;
    PhysiK::Tet tetA;
    tetA.node0 = 0;
    tetA.node1 = 1;
    tetA.node2 = 2;
    tetA.node3 = 3;
    PhysiK::Tet tetB;
    tetB.node0 = 1;
    tetB.node1 = 2;
    tetB.node2 = 3;
    tetB.node3 = 4;

    matrix.BuildPattern(5, BuildSparsePatternFromTetConnectivity({tetA, tetB}));

    assert(matrix.blockCount == 5);
    assert(matrix.values.size() == 23);
    assert(HasSparseBlock(matrix, 1, 2));
    assert(HasSparseBlock(matrix, 2, 1));
    assert(HasSparseBlock(matrix, 4, 4));
    assert(!HasSparseBlock(matrix, 0, 4));
}

void TetMeshComponentCachesFemSparsePattern()
{
    PhysiK::TetMeshPhysicsComponent component;
    component.localToGlobalNodeIndex = {0, 1, 2, 3};
    PhysiK::Tet tet = CreateUnitTet();
    component.tets.push_back(tet);
    component.EnsureFemSparsePattern(4);

    assert(!component.femSparsePatternDirty);
    assert(component.GetFemSparseMatrix().blockCount == 4);
    assert(component.GetFemSparseMatrix().values.size() == 16);

    component.EnsureFemSparsePattern(4);
    assert(component.GetFemSparseMatrix().values.size() == 16);

    component.MarkFemSparsePatternDirty();
    assert(component.femSparsePatternDirty);
    component.EnsureFemSparsePattern(5);
    assert(!component.femSparsePatternDirty);
    assert(component.GetFemSparseMatrix().blockCount == 5);
}

void TetMeshComponentMapsLocalFemPatternToGlobalSolverNodes()
{
    PhysiK::TetMeshPhysicsComponent component;
    component.localToGlobalNodeIndex = {2, 4, 5, 7};
    component.tets.push_back(CreateUnitTet());
    component.EnsureFemSparsePattern(8);

    const PhysiK::SparseBlockMatrix& matrix = component.GetFemSparseMatrix();
    assert(matrix.blockCount == 8);
    assert(matrix.values.size() == 16);
    assert(HasSparseBlock(matrix, 2, 4));
    assert(HasSparseBlock(matrix, 7, 5));
    assert(!HasSparseBlock(matrix, 0, 1));
    assert(!HasSparseBlock(matrix, 1, 1));
}

void ConjugateGradientSolvesDiagonalSparseSystem()
{
    PhysiK::SparseBlockMatrix matrix;
    matrix.BuildPattern(1, {{0, 0}});
    matrix.AddBlock(
        0,
        0,
        PhysiK::Mat3::FromColumns(
            PhysiK::Vec3{4.0f, 0.0f, 0.0f},
            PhysiK::Vec3{0.0f, 9.0f, 0.0f},
            PhysiK::Vec3{0.0f, 0.0f, 16.0f}));

    const std::vector<float> rhs = {4.0f, 18.0f, 48.0f};
    std::vector<float> solution;
    PhysiK::ConjugateGradientSettings settings;
    settings.maxIterations = 16;
    settings.tolerance = 1.0e-6f;
    settings.useJacobiPreconditioner = true;

    const PhysiK::ConjugateGradientResult result =
        PhysiK::SolveConjugateGradient(matrix, rhs, solution, settings);

    assert(result.converged);
    assert(solution.size() == 3);
    assert(NearlyEqual(solution[0], 1.0f, 0.0001f));
    assert(NearlyEqual(solution[1], 2.0f, 0.0001f));
    assert(NearlyEqual(solution[2], 3.0f, 0.0001f));
}

void ConjugateGradientSolvesCoupledSparseSystem()
{
    PhysiK::SparseBlockMatrix matrix;
    matrix.BuildPattern(2, {{0, 0}, {0, 1}, {1, 0}, {1, 1}});
    matrix.AddBlock(0, 0, DiagonalBlock(4.0f));
    matrix.AddBlock(0, 1, DiagonalBlock(-1.0f));
    matrix.AddBlock(1, 0, DiagonalBlock(-1.0f));
    matrix.AddBlock(1, 1, DiagonalBlock(3.0f));

    const std::vector<float> expected = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f};
    std::vector<float> rhs;
    matrix.Multiply(expected, rhs);

    std::vector<float> solution;
    PhysiK::ConjugateGradientSettings settings;
    settings.maxIterations = 32;
    settings.tolerance = 1.0e-6f;
    settings.useJacobiPreconditioner = true;

    const PhysiK::ConjugateGradientResult result =
        PhysiK::SolveConjugateGradient(matrix, rhs, solution, settings);

    assert(result.converged);
    assert(solution.size() == expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i)
    {
        assert(NearlyEqual(solution[i], expected[i], 0.0001f));
    }
}

void CurrentLinearSolverSolvesKnownSparseSystem()
{
    PhysiK::SparseBlockMatrix matrix;
    matrix.BuildPattern(1, {{0, 0}});
    matrix.AddBlock(
        0,
        0,
        PhysiK::Mat3::FromColumns(
            PhysiK::Vec3{2.0f, 0.0f, 0.0f},
            PhysiK::Vec3{0.0f, 3.0f, 0.0f},
            PhysiK::Vec3{0.0f, 0.0f, 4.0f}));

    const std::vector<float> rhs = {2.0f, 6.0f, 12.0f};
    std::vector<float> solution;

    PhysiK::CurrentLinearSolver solver;
    PhysiK::LinearSolveSettings settings;
    settings.maxIterations = 16;
    settings.tolerance = 1.0e-6f;
    settings.useJacobiPreconditioner = true;

    const PhysiK::LinearSolveResult result =
        solver.Solve(matrix, rhs, solution, settings);

    assert(result.converged);
    assert(solution.size() == rhs.size());
    assert(NearlyEqual(solution[0], 1.0f, 0.0001f));
    assert(NearlyEqual(solution[1], 2.0f, 0.0001f));
    assert(NearlyEqual(solution[2], 3.0f, 0.0001f));
}

void SolverDataFailedImplicitSolveLeavesNoDeltaVelocity()
{
    std::vector<PhysiK::Node> nodes(1);
    nodes[0].position = PhysiK::Vec3{0.0f, 0.0f, 0.0f};
    nodes[0].fixed = true;

    PhysiK::SolverData solverData;
    solverData.AddNodeMass(0, 1.0f);
    solverData.AddNodeForce(0, PhysiK::Vec3{1.0f, 0.0f, 0.0f});

    assert(!solverData.PrecomputeImplicitSolve(nodes, 0.01f));
    assert(!solverData.SolveImplicitLinearSystem());
    assert(solverData.GetDeltaVelocity().empty());
    assert(solverData.GetDynamicBlockForNode(0) < 0);
}

void PerformanceLoggingWritesCsvForImplicitStep()
{
#if defined(PHYSIK_ENABLE_PERF_LOGGING)
    const std::filesystem::path logPath =
        std::filesystem::path("logs") / "physik_performance_test.csv";
    std::error_code error;
    std::filesystem::remove(logPath, error);

    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    Require(world != nullptr, "world creation failed for performance logging test");

    int nodes[4] = {};
    CreateSingleTetWithMaterial(world, nodes, 100.0f, 0.0f, 1.0f);
    PHYSIK_SetSolverMode(world, 1);
    PHYSIK_SetPerformanceLogPath(world, logPath.string().c_str());

    PHYSIK_Step(world, 0.01f);
    PHYSIK_DestroyWorld(world);

    Require(
        std::filesystem::exists(logPath),
        "performance logging did not create the CSV file");

    std::ifstream file(logPath);
    std::string header;
    std::string row;
    std::getline(file, header);
    std::getline(file, row);

    Require(
        header.find("frameIndex,substepIndex,dt,totalStepMs") != std::string::npos,
        "performance CSV header is missing expected timing columns");
    Require(
        header.find("cgIterations,cgResidual,dynamicBlockCount,tetCount,activeTetCount") !=
            std::string::npos,
        "performance CSV header is missing expected solver/topology columns");
    Require(!row.empty(), "performance CSV did not contain a data row");

    std::filesystem::remove(logPath, error);
#else
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    Require(world != nullptr, "world creation failed for performance logging no-op test");

    PHYSIK_SetPerformanceLogPath(world, "logs/physik_performance_disabled_test.csv");
    PHYSIK_EnablePerformanceLogging(world, 1);
    PHYSIK_Step(world, 0.01f);

    PHYSIK_DestroyWorld(world);
#endif
}

void ImplicitEulerLinearTetUsesSparseCgPath()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    int nodes[4] = {};
    CreateSingleTetWithMaterial(world, nodes, 100.0f, 0.0f, 1.0f);
    PHYSIK_SetSolverMode(world, 1);
    PHYSIK_SetGravity(world, 0.0f, -9.81f, 0.0f);

    PHYSIK_Step(world, 0.01f);
    const Point velocity = GetNodeVelocity(world, nodes[3]);
    const Point position = GetNodePosition(world, nodes[3]);

    assert(IsFinite(velocity.x));
    assert(IsFinite(velocity.y));
    assert(IsFinite(velocity.z));
    assert(IsFinite(position.x));
    assert(IsFinite(position.y));
    assert(IsFinite(position.z));

    PHYSIK_DestroyWorld(world);
}

void ImplicitEulerCorotationalTetUsesSparseCgPath()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    const int node0 = AddNode(world, 0.0f, 0.0f, 0.0f);
    const int node1 = AddNode(world, 1.0f, 0.0f, 0.0f);
    const int node2 = AddNode(world, 0.0f, 1.0f, 0.0f);
    const int node3 = AddNode(world, 0.0f, 0.0f, 1.0f);
    const int nodes[] = {node0, node1, node2, node3};
    const int tetNodeIndices[] = {node0, node1, node2, node3};
    PhysikMaterialDesc material = MakeMaterialDesc(1.0f, 100.0f, 0.3f, 0.0f);
    const PhysiK::ComponentHandle tetMesh =
        PHYSIK_CreateTetMeshPhysicsComponent(world, nodes, 4, tetNodeIndices, 1, &material, 1);
    assert(PHYSIK_IsComponentHandleValid(world, tetMesh) == 1);

    PHYSIK_SetSolverMode(world, 1);
    PHYSIK_SetGravity(world, 0.0f, -9.81f, 0.0f);
    PHYSIK_Step(world, 0.01f);

    const Point velocity = GetNodeVelocity(world, node3);
    const Point position = GetNodePosition(world, node3);
    assert(IsFinite(velocity.x));
    assert(IsFinite(velocity.y));
    assert(IsFinite(velocity.z));
    assert(IsFinite(position.x));
    assert(IsFinite(position.y));
    assert(IsFinite(position.z));

    PHYSIK_DestroyWorld(world);
}

void MultiTetImplicitEulerSparseCgSmokeTest()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    const int nodes[] = {
        AddNode(world, 0.0f, 0.0f, 0.0f),
        AddNode(world, 1.0f, 0.0f, 0.0f),
        AddNode(world, 0.0f, 1.0f, 0.0f),
        AddNode(world, 0.0f, 0.0f, 1.0f),
        AddNode(world, 1.0f, 1.0f, 1.0f),
        AddNode(world, 1.0f, 0.0f, 1.0f),
        AddNode(world, 0.0f, 1.0f, 1.0f)};
    const int tetNodeIndices[] = {
        nodes[0], nodes[1], nodes[2], nodes[3],
        nodes[1], nodes[2], nodes[3], nodes[4],
        nodes[1], nodes[3], nodes[4], nodes[5],
        nodes[2], nodes[3], nodes[4], nodes[6]};
    PhysikMaterialDesc material = MakeMaterialDesc(1.0f, 75.0f, 0.3f, 0.0f);
    const PhysiK::ComponentHandle tetMesh =
        PHYSIK_CreateTetMeshPhysicsComponent(world, nodes, 7, tetNodeIndices, 4, &material, 0);
    assert(PHYSIK_IsComponentHandleValid(world, tetMesh) == 1);

    PHYSIK_SetSolverMode(world, 1);
    PHYSIK_SetGravity(world, 0.0f, -9.81f, 0.0f);
    PHYSIK_Step(world, 0.01f);

    for (int node : nodes)
    {
        const Point position = GetNodePosition(world, node);
        const Point velocity = GetNodeVelocity(world, node);
        assert(IsFinite(position.x));
        assert(IsFinite(position.y));
        assert(IsFinite(position.z));
        assert(IsFinite(velocity.x));
        assert(IsFinite(velocity.y));
        assert(IsFinite(velocity.z));
    }

    PHYSIK_DestroyWorld(world);
}

void TetActiveStateDefaultsAndNoOpsAreSafe()
{
    PhysiK::TetMeshComponent component;
    component.tets.resize(2);

    assert(component.IsTetActive(0));
    assert(component.IsTetActive(1));
    assert(!component.IsTetActive(-1));
    assert(!component.IsTetActive(2));
    assert(component.GetActiveTetCount() == 2);

    component.DeactivateTet(0);
    assert(!component.IsTetActive(0));
    assert(component.IsTetActive(1));
    assert(component.GetActiveTetCount() == 1);

    component.DeactivateTet(0);
    component.DeactivateTet(-1);
    component.DeactivateTet(4);
    assert(component.GetActiveTetCount() == 1);

    component.SetTetActive(0, true);
    assert(component.IsTetActive(0));
    assert(component.GetActiveTetCount() == 2);
}

void TetMeshComponentStoresGeometryWithoutWorldNodes()
{
    PhysiK::TetMeshComponent component;
    const PhysiK::Vec3 positions[] = {
        PhysiK::Vec3{0.0f, 0.0f, 0.0f},
        PhysiK::Vec3{1.0f, 0.0f, 0.0f},
        PhysiK::Vec3{0.0f, 1.0f, 0.0f},
        PhysiK::Vec3{0.0f, 0.0f, 1.0f}};
    const int tetIndices[] = {0, 1, 2, 3};

    const PhysiK::GeneratedTetMesh generatedMesh =
        PhysiK::TetMeshGenerator::Generate(positions, 4, tetIndices, 1);
    component.SetGeometry(generatedMesh);

    assert(component.GetNodeCount() == 4);
    assert(component.GetTetCount() == 1);
    assert(component.GetTetNodeIndex(0, 2) == 2);
    assert(component.GetGlobalNodeIndex(0) == -1);
    assert(NearlyEqual(component.GetLocalRestPosition(3), positions[3]));
    assert(NearlyEqual(component.GetLocalCurrentPosition(3), positions[3]));

    component.SetLocalCurrentPosition(3, PhysiK::Vec3{0.0f, 0.0f, 2.0f});
    assert(NearlyEqual(
        component.GetLocalCurrentPosition(3),
        PhysiK::Vec3{0.0f, 0.0f, 2.0f}));
}

void TetActiveStateIsExposedThroughNativeApi()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    int nodes[5] = {};
    const PhysiK::ComponentHandle tetMesh = CreateTwoTetMesh(world, nodes);

    assert(PHYSIK_GetTetMeshTetCount(world, tetMesh) == 2);
    assert(PHYSIK_GetActiveTetCount(world, tetMesh) == 2);
    assert(PHYSIK_IsTetActive(world, tetMesh, 0) == 1);
    assert(PHYSIK_IsTetActive(world, tetMesh, 1) == 1);
    assert(PHYSIK_IsTetActive(world, tetMesh, 2) == 0);

    PHYSIK_DeactivateTet(world, tetMesh, 1);
    assert(PHYSIK_IsTetActive(world, tetMesh, 0) == 1);
    assert(PHYSIK_IsTetActive(world, tetMesh, 1) == 0);
    assert(PHYSIK_GetActiveTetCount(world, tetMesh) == 1);

    PHYSIK_DeactivateTet(world, tetMesh, 1);
    PHYSIK_DeactivateTet(world, tetMesh, -1);
    PHYSIK_DeactivateTet(world, tetMesh, 3);
    assert(PHYSIK_GetActiveTetCount(world, tetMesh) == 1);

    PHYSIK_SetTetActive(world, tetMesh, 1, 1);
    assert(PHYSIK_IsTetActive(world, tetMesh, 1) == 1);
    assert(PHYSIK_GetActiveTetCount(world, tetMesh) == 2);

    PHYSIK_DestroyWorld(world);
}

void InactiveTetsAreSkippedByFemForceAndStiffnessAssembly()
{
    std::vector<PhysiK::Node> nodes = CreateTwoTetNodes();
    PhysiK::Tet activeTet = CreateUnitTet(200.0f);
    PhysiK::Tet inactiveTet = CreateLowerUnitTet(200.0f);
    PhysiK::FEMModel::InitializeTetRestData(activeTet, RestPositionsFromNodes(nodes));
    PhysiK::FEMModel::InitializeTetRestData(inactiveTet, RestPositionsFromNodes(nodes));
    inactiveTet.active = false;
    nodes[4].position.z -= 0.25f;

    for (PhysiK::FemModel model : {PhysiK::FemModel::Linear, PhysiK::FemModel::Corotational})
    {
        PhysiK::SolverData solverData;
        const bool implemented = ComputeFemForcesIntoSolverData(
            model,
            {activeTet, inactiveTet},
            nodes,
            solverData);

        assert(implemented);
        assert(LengthSquared(SumForcesForNode(solverData, 4)) < 0.000001f);
        assert(FindBlock(solverData, 4, 4) == nullptr);
        assert(FindBlock(solverData, 0, 4) == nullptr);
        assert(!solverData.GetStiffnessBlocks().empty());
    }
}

void DeactivatedTetsAreSkippedByLumpedMassAssembly()
{
    auto run = [](bool deactivateSecondTet)
    {
        PhysiK::WorldHandle world = PHYSIK_CreateWorld();
        assert(world != nullptr);

        int nodes[5] = {};
        const PhysiK::ComponentHandle tetMesh = CreateTwoTetMesh(world, nodes, 0.0f, 1.0f);
        PHYSIK_SetGravity(world, 0.0f, 0.0f, 0.0f);
        if (deactivateSecondTet)
        {
            PHYSIK_DeactivateTet(world, tetMesh, 1);
            assert(PHYSIK_GetActiveTetCount(world, tetMesh) == 1);
        }

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
            0.0f,
            0.0f,
            1.0f,
            1.0f,
            0.0f);
        PHYSIK_Step(world, 0.1f);
        const Point velocity = GetNodeVelocity(world, nodes[0]);
        PHYSIK_DestroyWorld(world);
        return velocity.z;
    };

    const float velocityWithBothTets = run(false);
    const float velocityWithOneTet = run(true);

    assert(velocityWithBothTets > 0.0f);
    assert(velocityWithOneTet > velocityWithBothTets * 1.5f);
}

void DeactivatingTetDoesNotDirtySparsePattern()
{
    PhysiK::TetMeshPhysicsComponent component;
    component.localToGlobalNodeIndex = {0, 1, 2, 3, 4};
    component.tets.push_back(CreateUnitTet());
    component.tets.push_back(CreateLowerUnitTet());
    component.tets.push_back(CreateUnitTet());
    component.tets.push_back(CreateLowerUnitTet());
    component.EnsureFemSparsePattern(5);

    const std::vector<int> rowStart = component.GetFemSparseMatrix().rowStart;
    const std::vector<int> colIndex = component.GetFemSparseMatrix().colIndex;
    const std::size_t valueCount = component.GetFemSparseMatrix().values.size();
    assert(!component.femSparsePatternDirty);

    component.DeactivateTet(1);
    assert(!component.femSparsePatternDirty);
    component.EnsureFemSparsePattern(5);

    assert(component.GetFemSparseMatrix().rowStart == rowStart);
    assert(component.GetFemSparseMatrix().colIndex == colIndex);
    assert(component.GetFemSparseMatrix().values.size() == valueCount);
}

void SmallTetMeshSimulatesAfterTetDeactivation()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    int nodes[5] = {};
    const PhysiK::ComponentHandle tetMesh = CreateTwoTetMesh(world, nodes, 25.0f, 1.0f);
    PHYSIK_SetSolverMode(world, 1);
    PHYSIK_SetGravity(world, 0.0f, -9.81f, 0.0f);
    PHYSIK_DeactivateTet(world, tetMesh, 1);

    PHYSIK_Step(world, 0.02f);
    const Point position = GetNodePosition(world, nodes[3]);
    const Point velocity = GetNodeVelocity(world, nodes[3]);

    assert(IsFinite(position.x));
    assert(IsFinite(position.y));
    assert(IsFinite(position.z));
    assert(IsFinite(velocity.x));
    assert(IsFinite(velocity.y));
    assert(IsFinite(velocity.z));
    assert(PHYSIK_GetActiveTetCount(world, tetMesh) == 1);

    PHYSIK_DestroyWorld(world);
}

void TetMeshComponentDefaultFemModelIsLinear()
{
    PhysiK::TetMeshPhysicsComponent component;

    assert(component.GetFemModel() == PhysiK::FemModel::Linear);
}

void TetMeshComponentStoresSelectedFemModel()
{
    PhysiK::TetMeshPhysicsComponent component;

    component.SetFemModel(PhysiK::FemModel::Linear);
    assert(component.GetFemModel() == PhysiK::FemModel::Linear);

    component.SetFemModel(PhysiK::FemModel::Corotational);
    assert(component.GetFemModel() == PhysiK::FemModel::Corotational);

    component.SetFemModel(PhysiK::FemModel::NeoHookean);
    assert(component.GetFemModel() == PhysiK::FemModel::NeoHookean);
}

void EventSystemDeliversSubscribedEventsOnlyOnce()
{
    PhysiK::EventSystem eventSystem;
    RecordingEventComponent listener;
    PhysiK::Component sender;
    listener.listenedEvents.push_back(PhysiK::PhysicsEventType::TetMeshTopologyChanged);
    sender.emittedEvents.push_back(PhysiK::PhysicsEventType::TetMeshTopologyChanged);

    eventSystem.Subscribe(&listener, PhysiK::PhysicsEventType::TetMeshTopologyChanged);
    eventSystem.Subscribe(&listener, PhysiK::PhysicsEventType::TetMeshTopologyChanged);

    const PhysiK::PhysicsEvent event{
        PhysiK::PhysicsEventType::TetMeshTopologyChanged,
        nullptr,
        &sender};
    eventSystem.Emit(event);

    assert(listener.eventCount == 1);
    assert(listener.lastEvent.type == PhysiK::PhysicsEventType::TetMeshTopologyChanged);
    assert(listener.lastEvent.sender == &sender);
    assert(listener.listenedEvents.size() == 1);
    assert(sender.emittedEvents.size() == 1);

    eventSystem.Unsubscribe(&listener, PhysiK::PhysicsEventType::TetMeshTopologyChanged);
    eventSystem.Emit(event);
    assert(listener.eventCount == 1);
}

void TopologyMeshComponentDeclaresEventsAndClearsDirtyFlag()
{
    PhysiK::WorldHandle worldHandle = PHYSIK_CreateWorld();
    assert(worldHandle != nullptr);
    PhysiK::World& world = *static_cast<PhysiK::World*>(worldHandle);

    int nodes[4] = {};
    const PhysiK::ComponentHandle tetMesh = CreateSingleTetMesh(worldHandle, nodes);

    PhysiK::TopologyMeshComponent topology(tetMesh);

    assert(!topology.listenedEvents.empty());
    assert(topology.listenedEvents.size() == 1);
    assert(topology.listenedEvents[0] ==
        PhysiK::PhysicsEventType::TetMeshTopologyChanged);
    assert(topology.emittedEvents.size() == 1);
    assert(topology.emittedEvents[0] ==
        PhysiK::PhysicsEventType::TopologyMeshUpdated);
    assert(topology.topologyDirty);

    topology.PostUpdate(world, 0.0f);
    assert(!topology.topologyDirty);
    assert(topology.GetIslandCount() == 1);

    PHYSIK_DestroyWorld(worldHandle);
}

void TopologyMeshComponentBuildsActiveTetIslands()
{
    PhysiK::WorldHandle worldHandle = PHYSIK_CreateWorld();
    assert(worldHandle != nullptr);
    PhysiK::World& world = *static_cast<PhysiK::World*>(worldHandle);

    const int node0 = AddNode(worldHandle, 0.0f, 0.0f, 0.0f);
    const int node1 = AddNode(worldHandle, 1.0f, 0.0f, 0.0f);
    const int node2 = AddNode(worldHandle, 0.0f, 1.0f, 0.0f);
    const int node3 = AddNode(worldHandle, 0.0f, 0.0f, 1.0f);
    const int node4 = AddNode(worldHandle, 0.0f, 0.0f, -1.0f);
    const int node5 = AddNode(worldHandle, 4.0f, 0.0f, 0.0f);
    const int node6 = AddNode(worldHandle, 5.0f, 0.0f, 0.0f);
    const int node7 = AddNode(worldHandle, 4.0f, 1.0f, 0.0f);
    const int node8 = AddNode(worldHandle, 4.0f, 0.0f, 1.0f);
    const int nodes[] = {
        node0, node1, node2, node3, node4, node5, node6, node7, node8};
    const int tetNodeIndices[] = {
        node0, node1, node2, node3,
        node0, node1, node2, node4,
        node5, node6, node7, node8};
    PhysikMaterialDesc material = MakeMaterialDesc(1.0f, 0.0f);

    const PhysiK::ComponentHandle tetMesh =
        PHYSIK_CreateTetMeshPhysicsComponent(
            worldHandle,
            nodes,
            9,
            tetNodeIndices,
            3,
            &material,
            0);
    assert(PHYSIK_IsComponentHandleValid(worldHandle, tetMesh) == 1);

    PhysiK::TopologyMeshComponent topology(tetMesh);
    topology.PostUpdate(world, 0.0f);

    assert(topology.GetIslandCount() == 2);
    assert(topology.GetTetIslandId(0) == 0);
    assert(topology.GetTetIslandId(1) == 0);
    assert(topology.GetTetIslandId(2) == 1);
    assert(topology.GetTetIslandId(-1) == -1);
    assert(topology.GetTetIslandId(3) == -1);

    PHYSIK_DeactivateTet(worldHandle, tetMesh, 0);
    topology.topologyDirty = true;
    topology.PostUpdate(world, 0.0f);

    assert(topology.GetIslandCount() == 2);
    assert(topology.GetTetIslandId(0) == -1);
    assert(topology.GetTetIslandId(1) == 0);
    assert(topology.GetTetIslandId(2) == 1);

    PHYSIK_DestroyWorld(worldHandle);
}

void SurfaceExtractionComponentExtractsActiveBoundaryFaces()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    int nodes[5] = {};
    const PhysiK::ComponentHandle tetMesh = CreateTwoTetMesh(world, nodes);
    PhysiK::SurfaceExtractionComponent surface(tetMesh);

    surface.RebuildSurface(*static_cast<PhysiK::World*>(world));
    assert(surface.GetSurfaceTriangleIndices().size() == 18);

    const std::vector<int>& indices = surface.GetSurfaceTriangleIndices();
    for (std::size_t index = 0; index < indices.size(); ++index)
    {
        assert(indices[index] >= 0);
    }

    PHYSIK_DeactivateTet(world, tetMesh, 0);
    surface.surfaceDirty = true;
    surface.PostUpdate(*static_cast<PhysiK::World*>(world), 0.0f);

    assert(!surface.surfaceDirty);
    assert(surface.GetSurfaceTriangleIndices().size() == 12);

    PHYSIK_DestroyWorld(world);
}

void SurfaceExtractionComponentWindsBoundaryFacesOutward()
{
    PhysiK::WorldHandle worldHandle = PHYSIK_CreateWorld();
    assert(worldHandle != nullptr);

    PhysiK::Vec3 positions[] = {
        PhysiK::Vec3{0.0f, 0.0f, 0.0f},
        PhysiK::Vec3{1.0f, 0.0f, 0.0f},
        PhysiK::Vec3{0.0f, 1.0f, 0.0f},
        PhysiK::Vec3{0.0f, 0.0f, 1.0f}};
    const int tetIndices[] = {0, 1, 2, 3};
    const PhysiK::ComponentHandle tetMeshHandle =
        PHYSIK_CreateTetMeshComponent(
            worldHandle,
            positions,
            4,
            tetIndices,
            1);
    assert(PHYSIK_IsComponentHandleValid(worldHandle, tetMeshHandle) == 1);

    const PhysiK::ComponentHandle surface =
        PHYSIK_CreateSurfaceExtractionComponent(worldHandle, tetMeshHandle);
    assert(PHYSIK_IsComponentHandleValid(worldHandle, surface) == 1);
    PHYSIK_Step(worldHandle, 0.0f);

    int indices[12] = {};
    assert(PHYSIK_CopySurfaceTriangleIndices(
        worldHandle,
        surface,
        indices,
        12) == 12);

    for (int index = 0; index < 12; index += 3)
    {
        const int node0 = indices[index + 0u];
        const int node1 = indices[index + 1u];
        const int node2 = indices[index + 2u];
        int oppositeNode = -1;
        for (int candidate = 0; candidate < 4; ++candidate)
        {
            if (candidate != node0 &&
                candidate != node1 &&
                candidate != node2)
            {
                oppositeNode = candidate;
                break;
            }
        }

        assert(oppositeNode >= 0);
        assert(TriangleWindsAwayFromOppositePoint(
            positions,
            node0,
            node1,
            node2,
            oppositeNode));
    }

    PHYSIK_DestroyWorld(worldHandle);
}

void TetMeshGeneratorWeldsNodesAndDropsDegenerateTets()
{
    const PhysiK::Vec3 positions[] = {
        PhysiK::Vec3{0.0f, 0.0f, 0.0f},
        PhysiK::Vec3{0.0f, 0.0f, 0.0f},
        PhysiK::Vec3{1.0f, 0.0f, 0.0f},
        PhysiK::Vec3{0.0f, 1.0f, 0.0f},
        PhysiK::Vec3{0.0f, 0.0f, 1.0f}};
    const int tetIndices[] = {
        0, 1, 2, 3,
        0, 2, 3, 4};

    const PhysiK::GeneratedTetMesh result =
        PhysiK::TetMeshGenerator::Generate(
            positions,
            5,
            tetIndices,
            2);

    assert(result.rawNodeCount == 5);
    assert(result.weldedNodeCount == 4);
    assert(result.weldedAwayNodeCount == 1);
    assert(result.rawTetCount == 2);
    assert(result.removedDegenerateTetCount == 1);
    assert(result.tetLocalNodeIndices.size() == 4);
    assert(result.topologyDiagnostics.boundaryFaceCount == 4);
    assert(result.topologyDiagnostics.internalFaceCount == 0);
    assert(result.topologyDiagnostics.nonManifoldFaceCount == 0);

}

void GeneratedTetMeshApiStoresCleanGeometry()
{
    const PhysiK::Vec3 positions[] = {
        PhysiK::Vec3{0.0f, 0.0f, 0.0f},
        PhysiK::Vec3{1.0f, 0.0f, 0.0f},
        PhysiK::Vec3{0.0f, 1.0f, 0.0f},
        PhysiK::Vec3{0.0f, 0.0f, 1.0f}};
    const int tetIndices[] = {0, 1, 2, 3};

    const PhysiK::GeneratedTetMeshHandle generatedMesh =
        PHYSIK_GenerateTetMesh(positions, 4, tetIndices, 1);
    assert(PHYSIK_IsGeneratedTetMeshHandleValid(generatedMesh) == 1);
    assert(PHYSIK_GetGeneratedTetMeshVertexCount(generatedMesh) == 4);
    assert(PHYSIK_GetGeneratedTetMeshTetCount(generatedMesh) == 1);
    assert(PHYSIK_GetGeneratedTetMeshTetIndexCount(generatedMesh) == 4);

    float x = -1.0f;
    float y = -1.0f;
    float z = -1.0f;
    assert(PHYSIK_GetGeneratedTetMeshVertex(
        generatedMesh,
        3,
        &x,
        &y,
        &z) == 1);
    assert(NearlyEqual(x, 0.0f));
    assert(NearlyEqual(y, 0.0f));
    assert(NearlyEqual(z, 1.0f));

    int nodeIndex = -1;
    assert(PHYSIK_GetGeneratedTetMeshTetNodeIndex(
        generatedMesh,
        2,
        &nodeIndex) == 1);
    assert(nodeIndex == 2);
    assert(PHYSIK_GetGeneratedTetMeshTetNodeIndex(
        generatedMesh,
        4,
        &nodeIndex) == 0);
    assert(PHYSIK_GetGeneratedTetMeshVertex(
        generatedMesh,
        -1,
        &x,
        &y,
        &z) == 0);

    PHYSIK_DestroyGeneratedTetMesh(generatedMesh);
    assert(PHYSIK_IsGeneratedTetMeshHandleValid(generatedMesh) == 0);
}

void TetMeshComponentCanBeCreatedFromGeneratedTetMeshApi()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    const PhysiK::Vec3 positions[] = {
        PhysiK::Vec3{0.0f, 0.0f, 0.0f},
        PhysiK::Vec3{1.0f, 0.0f, 0.0f},
        PhysiK::Vec3{0.0f, 1.0f, 0.0f},
        PhysiK::Vec3{0.0f, 0.0f, 1.0f}};
    const int tetIndices[] = {0, 1, 2, 3};
    const PhysiK::GeneratedTetMeshHandle generatedMesh =
        PHYSIK_GenerateTetMesh(positions, 4, tetIndices, 1);

    const PhysiK::ComponentHandle tetMesh =
        PHYSIK_CreateTetMeshComponentFromGeneratedTetMesh(
            world,
            generatedMesh);
    assert(PHYSIK_IsComponentHandleValid(world, tetMesh) == 1);
    assert(PHYSIK_GetTetMeshNodeCount(world, tetMesh) == 4);
    assert(PHYSIK_GetTetMeshTetCount(world, tetMesh) == 1);

    PHYSIK_DestroyGeneratedTetMesh(generatedMesh);
    PHYSIK_DestroyWorld(world);
}

void TetMeshPhysicsComponentCanBeCreatedFromGeneratedTetMeshApi()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    const PhysiK::Vec3 positions[] = {
        PhysiK::Vec3{0.0f, 0.0f, 0.0f},
        PhysiK::Vec3{1.0f, 0.0f, 0.0f},
        PhysiK::Vec3{0.0f, 1.0f, 0.0f},
        PhysiK::Vec3{0.0f, 0.0f, 1.0f}};
    const int tetIndices[] = {0, 1, 2, 3};
    const PhysiK::GeneratedTetMeshHandle generatedMesh =
        PHYSIK_GenerateTetMesh(positions, 4, tetIndices, 1);
    PhysikMaterialDesc material = MakeMaterialDesc(1.0f, 25.0f, 0.3f, 0.0f);

    const PhysiK::ComponentHandle tetMesh =
        PHYSIK_CreateTetMeshPhysicsComponentFromGeneratedTetMesh(
            world,
            generatedMesh,
            &material);
    assert(PHYSIK_IsComponentHandleValid(world, tetMesh) == 1);
    assert(PHYSIK_GetTetMeshNodeCount(world, tetMesh) == 4);
    assert(PHYSIK_GetTetMeshTetCount(world, tetMesh) == 1);
    assert(NearlyEqual(GetNodePosition(world, 3).z, 1.0f));

    PHYSIK_DestroyGeneratedTetMesh(generatedMesh);
    PHYSIK_DestroyWorld(world);
}

void GeneratedTetMeshApiWeldsDuplicateNodes()
{
    const PhysiK::Vec3 positions[] = {
        PhysiK::Vec3{0.0f, 0.0f, 0.0f},
        PhysiK::Vec3{0.0f, 0.0f, 0.0f},
        PhysiK::Vec3{1.0f, 0.0f, 0.0f},
        PhysiK::Vec3{0.0f, 1.0f, 0.0f},
        PhysiK::Vec3{0.0f, 0.0f, 1.0f}};
    const int tetIndices[] = {
        0, 1, 2, 3,
        0, 2, 3, 4};

    const PhysiK::GeneratedTetMeshHandle generatedMesh =
        PHYSIK_GenerateTetMesh(positions, 5, tetIndices, 2);
    assert(PHYSIK_IsGeneratedTetMeshHandleValid(generatedMesh) == 1);
    assert(PHYSIK_GetGeneratedTetMeshVertexCount(generatedMesh) == 4);
    assert(PHYSIK_GetGeneratedTetMeshTetCount(generatedMesh) == 1);
    assert(PHYSIK_GetGeneratedTetMeshTetIndexCount(generatedMesh) == 4);

    for (int i = 0; i < PHYSIK_GetGeneratedTetMeshTetIndexCount(generatedMesh); ++i)
    {
        int nodeIndex = -1;
        assert(PHYSIK_GetGeneratedTetMeshTetNodeIndex(
            generatedMesh,
            i,
            &nodeIndex) == 1);
        assert(nodeIndex >= 0);
        assert(nodeIndex < PHYSIK_GetGeneratedTetMeshVertexCount(generatedMesh));
    }

    PHYSIK_DestroyGeneratedTetMesh(generatedMesh);
}

void TetMeshPhysicsCanConsumeGeneratedTetMeshDirectly()
{
    PhysiK::WorldHandle worldHandle = PHYSIK_CreateWorld();
    assert(worldHandle != nullptr);
    PhysiK::World& world = *static_cast<PhysiK::World*>(worldHandle);

    PhysiK::GeneratedTetMesh generatedMesh;
    generatedMesh.positions = {
        PhysiK::Vec3{0.0f, 0.0f, 0.0f},
        PhysiK::Vec3{1.0f, 0.0f, 0.0f},
        PhysiK::Vec3{0.0f, 1.0f, 0.0f},
        PhysiK::Vec3{0.0f, 0.0f, 1.0f}};
    generatedMesh.tetLocalNodeIndices = {0, 1, 2, 3};

    PhysiK::TetMeshPhysicsComponentDesc desc;
    desc.material = PhysiK::Material{1.0f, 25.0f, 0.3f, 0.0f};
    desc.femModel = PhysiK::FemModel::Linear;

    std::unique_ptr<PhysiK::TetMeshPhysicsComponent> component =
        PhysiK::TetMeshPhysicsComponent::CreateFromGeneratedTetMesh(
            world,
            generatedMesh,
            desc);

    assert(component != nullptr);
    assert(component->GetNodeCount() == 4);
    assert(component->GetTetCount() == 1);
    assert(component->localToGlobalNodeIndex.size() == 4);
    assert(component->GetGlobalNodeIndex(0) == 0);
    assert(component->GetGlobalNodeIndex(3) == 3);
    assert(GetNodePosition(worldHandle, 3).z == 1.0f);

    PHYSIK_DestroyWorld(worldHandle);
}

void TetMeshCreationWeldsDuplicateSharedFaceNodes()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    const PhysiK::Vec3 positions[] = {
        PhysiK::Vec3{0.0f, 0.0f, 0.0f},
        PhysiK::Vec3{1.0f, 0.0f, 0.0f},
        PhysiK::Vec3{0.0f, 1.0f, 0.0f},
        PhysiK::Vec3{0.0f, 0.0f, 1.0f},
        PhysiK::Vec3{0.0f, 0.0f, 0.0f},
        PhysiK::Vec3{1.0f, 0.0f, 0.0f},
        PhysiK::Vec3{0.0f, 1.0f, 0.0f},
        PhysiK::Vec3{0.0f, 0.0f, -1.0f}};
    const int tetIndices[] = {
        0, 1, 2, 3,
        4, 5, 6, 7};

    const PhysiK::ComponentHandle tetMesh =
        PHYSIK_CreateTetMeshComponent(
            world,
            positions,
            8,
            tetIndices,
            2);
    assert(PHYSIK_IsComponentHandleValid(world, tetMesh) == 1);
    assert(PHYSIK_GetTetMeshNodeCount(world, tetMesh) == 5);
    assert(PHYSIK_GetTetMeshTetCount(world, tetMesh) == 2);

    const PhysiK::ComponentHandle surface =
        PHYSIK_CreateSurfaceExtractionComponent(world, tetMesh);
    assert(PHYSIK_IsComponentHandleValid(world, surface) == 1);
    PHYSIK_Step(world, 0.0f);

    assert(PHYSIK_GetSurfaceTriangleIndexCount(world, surface) == 18);

    PHYSIK_DestroyWorld(world);
}

void SurfaceExtractionComponentCanBeCreatedThroughNativeApi()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    const PhysiK::Vec3 positions[] = {
        PhysiK::Vec3{0.0f, 0.0f, 0.0f},
        PhysiK::Vec3{1.0f, 0.0f, 0.0f},
        PhysiK::Vec3{0.0f, 1.0f, 0.0f},
        PhysiK::Vec3{0.0f, 0.0f, 1.0f}};
    const int tetIndices[] = {0, 1, 2, 3};
    const PhysiK::ComponentHandle tetMesh =
        PHYSIK_CreateTetMeshComponent(
            world,
            positions,
            4,
            tetIndices,
            1);
    assert(PHYSIK_IsComponentHandleValid(world, tetMesh) == 1);

    const PhysiK::ComponentHandle surface =
        PHYSIK_CreateSurfaceExtractionComponent(world, tetMesh);

    assert(PHYSIK_IsComponentHandleValid(world, surface) == 1);
    assert(PHYSIK_GetSurfaceTriangleIndexCount(world, surface) == 0);

    PHYSIK_Step(world, 0.0f);
    assert(PHYSIK_GetSurfaceTriangleIndexCount(world, surface) == 12);

    int copiedIndices[12] = {};
    assert(PHYSIK_CopySurfaceTriangleIndices(
        world,
        surface,
        copiedIndices,
        12) == 12);
    for (int copiedIndex : copiedIndices)
    {
        assert(copiedIndex >= 0);
        assert(copiedIndex < 4);
    }

    assert(PHYSIK_GetSurfaceTriangleIndexCount(nullptr, surface) == 0);
    assert(PHYSIK_GetSurfaceTriangleIndexCount(world, tetMesh) == 0);
    assert(PHYSIK_CopySurfaceTriangleIndices(nullptr, surface, copiedIndices, 12) == 0);
    assert(PHYSIK_CopySurfaceTriangleIndices(world, tetMesh, copiedIndices, 12) == 0);
    assert(PHYSIK_CopySurfaceTriangleIndices(world, surface, nullptr, 12) == 0);
    assert(PHYSIK_CopySurfaceTriangleIndices(world, surface, copiedIndices, 0) == 0);

    assert(PHYSIK_CreateSurfaceExtractionComponent(nullptr, tetMesh).IsValid() == false);
    assert(PHYSIK_CreateSurfaceExtractionComponent(
        world,
        PhysiK::ComponentHandle{}).IsValid() == false);

    const PhysiK::ComponentHandle visual =
        PHYSIK_CreateVisualMeshComponent(world, tetMesh);
    assert(PHYSIK_IsComponentHandleValid(world, visual) == 1);
    assert(PHYSIK_CreateSurfaceExtractionComponent(world, visual).IsValid() == false);

    PHYSIK_DestroyWorld(world);
}

void SurfaceVisualComponentBuildsRenderReadySurface()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    const PhysiK::Vec3 positions[] = {
        PhysiK::Vec3{0.0f, 0.0f, 0.0f},
        PhysiK::Vec3{1.0f, 0.0f, 0.0f},
        PhysiK::Vec3{0.0f, 1.0f, 0.0f},
        PhysiK::Vec3{0.0f, 0.0f, 1.0f}};
    const int tetIndices[] = {0, 1, 2, 3};
    const PhysiK::ComponentHandle tetMesh =
        PHYSIK_CreateTetMeshComponent(
            world,
            positions,
            4,
            tetIndices,
            1);
    assert(PHYSIK_IsComponentHandleValid(world, tetMesh) == 1);

    const PhysiK::ComponentHandle surface =
        PHYSIK_CreateSurfaceExtractionComponent(world, tetMesh);
    assert(PHYSIK_IsComponentHandleValid(world, surface) == 1);

    const PhysiK::ComponentHandle visual =
        PHYSIK_CreateSurfaceVisualComponent(world, surface);
    assert(PHYSIK_IsComponentHandleValid(world, visual) == 1);
    assert(PHYSIK_CreateSurfaceVisualComponent(nullptr, surface).IsValid() == false);
    assert(PHYSIK_CreateSurfaceVisualComponent(world, tetMesh).IsValid() == false);

    assert(PHYSIK_SetTetMeshLocalCurrentPosition(
        world,
        tetMesh,
        0,
        0.0f,
        0.0f,
        2.0f) == 1);

    PHYSIK_Step(world, 0.0f);

    assert(PHYSIK_GetSurfaceVisualVertexCount(world, visual) == 12);
    assert(PHYSIK_GetSurfaceVisualTriangleIndexCount(world, visual) == 12);
    assert(PHYSIK_GetSurfaceVisualNormalCount(world, visual) == 12);
    assert(PHYSIK_GetSurfaceVisualVertexCount(nullptr, visual) == 0);
    assert(PHYSIK_GetSurfaceVisualTriangleIndexCount(world, tetMesh) == 0);
    assert(PHYSIK_GetSurfaceVisualNormalCount(world, surface) == 0);

    Point vertex;
    assert(PHYSIK_GetSurfaceVisualVertex(
        world,
        visual,
        0,
        &vertex.x,
        &vertex.y,
        &vertex.z) == 1);

    Point normal;
    assert(PHYSIK_GetSurfaceVisualNormal(
        world,
        visual,
        0,
        &normal.x,
        &normal.y,
        &normal.z) == 1);
    assert(NearlyEqual(
        std::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z),
        1.0f));

    int triangleIndex = -1;
    assert(PHYSIK_GetSurfaceVisualTriangleIndex(
        world,
        visual,
        0,
        &triangleIndex) == 1);
    assert(triangleIndex == 0);

    assert(PHYSIK_GetSurfaceVisualVertex(
        world,
        visual,
        -1,
        &vertex.x,
        &vertex.y,
        &vertex.z) == 0);
    assert(PHYSIK_GetSurfaceVisualVertex(
        world,
        visual,
        12,
        &vertex.x,
        &vertex.y,
        &vertex.z) == 0);
    assert(PHYSIK_GetSurfaceVisualVertex(
        world,
        visual,
        0,
        nullptr,
        &vertex.y,
        &vertex.z) == 0);
    assert(PHYSIK_GetSurfaceVisualNormal(
        world,
        visual,
        0,
        &normal.x,
        nullptr,
        &normal.z) == 0);
    assert(PHYSIK_GetSurfaceVisualTriangleIndex(
        world,
        visual,
        0,
        nullptr) == 0);

    bool foundMovedNode = false;
    for (int vertexIndex = 0; vertexIndex < PHYSIK_GetSurfaceVisualVertexCount(world, visual);
         ++vertexIndex)
    {
        assert(PHYSIK_GetSurfaceVisualVertex(
            world,
            visual,
            vertexIndex,
            &vertex.x,
            &vertex.y,
            &vertex.z) == 1);
        if (NearlyEqual(vertex.x, 0.0f) &&
            NearlyEqual(vertex.y, 0.0f) &&
            NearlyEqual(vertex.z, 2.0f))
        {
            foundMovedNode = true;
        }
    }

    assert(foundMovedNode);

    const Point local0 = GetTetMeshLocalCurrentPosition(world, tetMesh, 0);
    assert(NearlyEqual(local0.z, 2.0f));

    assert(PHYSIK_SetTetMeshLocalCurrentPosition(
        world,
        tetMesh,
        1,
        3.0f,
        0.0f,
        0.0f) == 1);
    PHYSIK_Step(world, 0.0f);
    assert(PHYSIK_GetSurfaceVisualVertexCount(world, visual) == 12);
    assert(PHYSIK_GetSurfaceVisualTriangleIndexCount(world, visual) == 12);
    assert(PHYSIK_GetSurfaceVisualNormalCount(world, visual) == 12);

    bool foundUpdatedNode = false;
    for (int vertexIndex = 0; vertexIndex < PHYSIK_GetSurfaceVisualVertexCount(world, visual);
         ++vertexIndex)
    {
        assert(PHYSIK_GetSurfaceVisualVertex(
            world,
            visual,
            vertexIndex,
            &vertex.x,
            &vertex.y,
            &vertex.z) == 1);
        if (NearlyEqual(vertex.x, 3.0f) &&
            NearlyEqual(vertex.y, 0.0f) &&
            NearlyEqual(vertex.z, 0.0f))
        {
            foundUpdatedNode = true;
        }
    }

    assert(foundUpdatedNode);

    PHYSIK_DestroyWorld(world);
}

void VisualMeshComponentDeclaresTopologyListenerAndClearsDirtyFlag()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    PhysiK::VisualMeshComponent visual;

    assert(visual.listenedEvents.size() == 1);
    assert(visual.listenedEvents[0] == PhysiK::PhysicsEventType::TetMeshTopologyChanged);

    visual.topologyDirty = true;
    visual.PostUpdate(*static_cast<PhysiK::World*>(world), 0.0f);
    assert(!visual.topologyDirty);

    PHYSIK_DestroyWorld(world);
}

void VisualMeshComponentCanBeCreatedThroughNativeApi()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    int nodes[4] = {};
    const PhysiK::ComponentHandle tetMesh = CreateSingleTetMesh(world, nodes);
    const PhysiK::ComponentHandle visual =
        PHYSIK_CreateVisualMeshComponent(world, tetMesh);

    assert(PHYSIK_IsComponentHandleValid(world, visual) == 1);

    PHYSIK_DeactivateTet(world, tetMesh, 0);
    PHYSIK_Step(world, 0.01f);

    PHYSIK_DestroyWorld(world);
}

void VisualMeshComponentStoresVisualMeshData()
{
    PhysiK::VisualMeshComponent visual;
    const PhysiK::Vec3 vertices[] = {
        PhysiK::Vec3{0.0f, 0.0f, 0.0f},
        PhysiK::Vec3{1.0f, 0.0f, 0.0f},
        PhysiK::Vec3{0.0f, 1.0f, 0.0f},
        PhysiK::Vec3{0.0f, 0.0f, 1.0f}};
    const int indices[] = {0, 1, 2, 0, 2, 3};

    visual.SetVisualMesh(vertices, 4, indices, 6);

    assert(visual.restVisualVertices.size() == 4);
    assert(visual.GetDeformedVertices().size() == 4);
    assert(visual.GetTriangleIndices().size() == 6);
    assert(visual.embeddedVertices.size() == 4);
    assert(visual.triangleValid.size() == 2);

    for (int i = 0; i < 4; ++i)
    {
        assert(NearlyEqual(visual.restVisualVertices[i], vertices[i]));
        assert(NearlyEqual(visual.GetDeformedVertices()[i], vertices[i]));
        assert(visual.embeddedVertices[i].tetIndex == -1);
        assert(!visual.embeddedVertices[i].valid);
    }

    for (int i = 0; i < 6; ++i)
    {
        assert(visual.GetTriangleIndices()[i] == indices[i]);
    }

    for (bool valid : visual.triangleValid)
    {
        assert(valid);
    }
}

void VisualMeshComponentBuildsBruteForceEmbedding()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    int nodes[4] = {};
    const PhysiK::ComponentHandle tetMesh = CreateSingleTetMesh(world, nodes);
    PhysiK::VisualMeshComponent visual(tetMesh, "test visual");
    const PhysiK::Vec3 vertices[] = {
        PhysiK::Vec3{0.25f, 0.25f, 0.25f},
        PhysiK::Vec3{2.0f, 2.0f, 2.0f}};

    visual.SetVisualMesh(vertices, 2, nullptr, 0);
    visual.BuildEmbedding(*static_cast<PhysiK::World*>(world));

    assert(visual.embeddedVertices.size() == 2);
    assert(visual.embeddedVertices[0].valid);
    assert(visual.embeddedVertices[0].tetIndex == 0);
    assert(NearlyEqual(visual.embeddedVertices[0].barycentric.x, 0.25f));
    assert(NearlyEqual(visual.embeddedVertices[0].barycentric.y, 0.25f));
    assert(NearlyEqual(visual.embeddedVertices[0].barycentric.z, 0.25f));
    assert(NearlyEqual(visual.embeddedVertices[0].barycentric.w, 0.25f));
    assert(!visual.embeddedVertices[1].valid);
    assert(visual.embeddedVertices[1].tetIndex == -1);

    PHYSIK_DeactivateTet(world, tetMesh, 0);
    visual.BuildEmbedding(*static_cast<PhysiK::World*>(world));
    assert(!visual.embeddedVertices[0].valid);
    assert(visual.embeddedVertices[0].tetIndex == -1);

    PHYSIK_DestroyWorld(world);
}

void VisualMeshComponentUpdatesDeformedVerticesFromHostTet()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    int nodes[4] = {};
    const PhysiK::ComponentHandle tetMesh = CreateSingleTetMesh(world, nodes);
    PhysiK::VisualMeshComponent visual(tetMesh, "test visual");
    const PhysiK::Vec3 vertices[] = {
        PhysiK::Vec3{0.25f, 0.25f, 0.25f},
        PhysiK::Vec3{2.0f, 2.0f, 2.0f}};

    visual.SetVisualMesh(vertices, 2, nullptr, 0);
    visual.BuildEmbedding(*static_cast<PhysiK::World*>(world));

    PHYSIK_SetNodePosition(world, nodes[0], 0.0f, 0.0f, 1.0f);
    PHYSIK_SetNodePosition(world, nodes[1], 2.0f, 0.0f, 1.0f);
    PHYSIK_SetNodePosition(world, nodes[2], 0.0f, 2.0f, 1.0f);
    PHYSIK_SetNodePosition(world, nodes[3], 0.0f, 0.0f, 3.0f);

    visual.UpdateDeformedVertices(*static_cast<PhysiK::World*>(world));

    assert(NearlyEqual(
        visual.GetDeformedVertices()[0],
        PhysiK::Vec3{0.5f, 0.5f, 1.5f}));
    assert(NearlyEqual(visual.GetDeformedVertices()[1], vertices[1]));

    visual.topologyDirty = true;
    visual.PostUpdate(*static_cast<PhysiK::World*>(world), 0.0f);
    assert(!visual.topologyDirty);
    assert(NearlyEqual(
        visual.GetDeformedVertices()[0],
        PhysiK::Vec3{0.5f, 0.5f, 1.5f}));

    PHYSIK_DeactivateTet(world, tetMesh, 0);
    PHYSIK_SetNodePosition(world, nodes[0], 10.0f, 10.0f, 10.0f);
    visual.UpdateDeformedVertices(*static_cast<PhysiK::World*>(world));
    assert(NearlyEqual(
        visual.GetDeformedVertices()[0],
        PhysiK::Vec3{0.5f, 0.5f, 1.5f}));

    PHYSIK_DestroyWorld(world);
}

void VisualMeshComponentCAPIExportsMeshBuffers()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    int nodes[4] = {};
    const PhysiK::ComponentHandle tetMesh = CreateSingleTetMesh(world, nodes);
    const PhysiK::ComponentHandle visualHandle =
        PHYSIK_CreateVisualMeshComponent(world, tetMesh);

    assert(PHYSIK_IsComponentHandleValid(world, visualHandle) == 1);
    assert(PHYSIK_GetVisualMeshVertexCount(world, visualHandle) == 0);
    assert(PHYSIK_GetVisualMeshTriangleIndexCount(world, visualHandle) == 0);

    const PhysiK::Vec3 vertices[] = {
        PhysiK::Vec3{0.25f, 0.25f, 0.25f},
        PhysiK::Vec3{0.5f, 0.25f, 0.25f},
        PhysiK::Vec3{0.25f, 0.5f, 0.25f},
        PhysiK::Vec3{2.0f, 2.0f, 2.0f}};
    const int indices[] = {0, 1, 2, 0, 2, 3};
    PHYSIK_SetVisualMeshData(world, visualHandle, vertices, 4, indices, 6);

    assert(PHYSIK_GetVisualMeshVertexCount(world, visualHandle) == 4);
    assert(PHYSIK_GetVisualMeshTriangleIndexCount(world, visualHandle) == 6);
    assert(PHYSIK_BuildVisualMeshEmbedding(world, visualHandle) == 1);
    assert(PHYSIK_GetVisualMeshTriangleIndexCount(world, visualHandle) == 3);

    PHYSIK_SetNodePosition(world, nodes[0], 0.0f, 0.0f, 1.0f);
    PHYSIK_SetNodePosition(world, nodes[1], 2.0f, 0.0f, 1.0f);
    PHYSIK_SetNodePosition(world, nodes[2], 0.0f, 2.0f, 1.0f);
    PHYSIK_SetNodePosition(world, nodes[3], 0.0f, 0.0f, 3.0f);
    PHYSIK_Step(world, 0.01f);

    PhysiK::Vec3 copiedVertices[4] = {};
    int copiedIndices[6] = {};
    assert(PHYSIK_CopyVisualMeshVertices(world, visualHandle, copiedVertices, 4) == 4);
    assert(NearlyEqual(copiedVertices[0], PhysiK::Vec3{0.5f, 0.5f, 1.5f}));
    assert(NearlyEqual(copiedVertices[1], PhysiK::Vec3{1.0f, 0.5f, 1.5f}, 0.1));
    assert(NearlyEqual(copiedVertices[2], PhysiK::Vec3{0.5f, 1.0f, 1.5f}, 0.1));
    assert(NearlyEqual(copiedVertices[3], vertices[3]));
    assert(PHYSIK_CopyVisualMeshTriangleIndices(world, visualHandle, copiedIndices, 6) == 3);
    assert(copiedIndices[0] == 0);
    assert(copiedIndices[1] == 1);
    assert(copiedIndices[2] == 2);

    PHYSIK_DeactivateTet(world, tetMesh, 0);
    PHYSIK_Step(world, 0.01f);
    assert(PHYSIK_GetVisualMeshTriangleIndexCount(world, visualHandle) == 0);
    assert(PHYSIK_CopyVisualMeshTriangleIndices(world, visualHandle, copiedIndices, 6) == 0);

    assert(PHYSIK_GetVisualMeshVertexCount(nullptr, visualHandle) == 0);
    assert(PHYSIK_GetVisualMeshTriangleIndexCount(world, tetMesh) == 0);
    assert(PHYSIK_CreateVisualMeshComponent(nullptr, tetMesh).IsValid() == false);
    assert(PHYSIK_CreateVisualMeshComponent(world, PhysiK::ComponentHandle{}).IsValid() == false);
    assert(PHYSIK_CreateVisualMeshComponent(world, visualHandle).IsValid() == false);
    assert(PHYSIK_BuildVisualMeshEmbedding(nullptr, visualHandle) == 0);
    assert(PHYSIK_BuildVisualMeshEmbedding(world, tetMesh) == 0);
    PHYSIK_SetVisualMeshData(nullptr, visualHandle, vertices, 4, indices, 6);
    PHYSIK_SetVisualMeshData(world, tetMesh, vertices, 4, indices, 6);
    assert(PHYSIK_CopyVisualMeshVertices(nullptr, visualHandle, copiedVertices, 3) == 0);
    assert(PHYSIK_CopyVisualMeshVertices(world, tetMesh, copiedVertices, 3) == 0);
    assert(PHYSIK_CopyVisualMeshVertices(world, visualHandle, nullptr, 3) == 0);
    assert(PHYSIK_CopyVisualMeshVertices(world, visualHandle, copiedVertices, 0) == 0);
    assert(PHYSIK_CopyVisualMeshTriangleIndices(world, visualHandle, nullptr, 3) == 0);
    assert(PHYSIK_CopyVisualMeshTriangleIndices(world, visualHandle, copiedIndices, 0) == 0);

    PHYSIK_DestroyWorld(world);
}

void TetMeshMapperComponentEmbedsDestinationAndFollowsSource()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    const PhysiK::Vec3 sourcePositions[] = {
        PhysiK::Vec3{0.0f, 0.0f, 0.0f},
        PhysiK::Vec3{1.0f, 0.0f, 0.0f},
        PhysiK::Vec3{0.0f, 1.0f, 0.0f},
        PhysiK::Vec3{0.0f, 0.0f, 1.0f}};
    const int sourceTetIndices[] = {0, 1, 2, 3};
    const PhysiK::ComponentHandle sourceTetMesh =
        PHYSIK_CreateTetMeshComponent(
            world,
            sourcePositions,
            4,
            sourceTetIndices,
            1);
    assert(PHYSIK_IsComponentHandleValid(world, sourceTetMesh) == 1);

    const PhysiK::Vec3 destinationPositions[] = {
        PhysiK::Vec3{0.25f, 0.25f, 0.25f},
        PhysiK::Vec3{0.50f, 0.25f, 0.25f},
        PhysiK::Vec3{0.25f, 0.50f, 0.25f},
        PhysiK::Vec3{0.25f, 0.25f, 0.50f},
        PhysiK::Vec3{2.0f, 2.0f, 2.0f}};
    const int destinationTetIndices[] = {0, 1, 2, 3};
    const PhysiK::ComponentHandle destinationTetMesh =
        PHYSIK_CreateTetMeshComponent(
            world,
            destinationPositions,
            5,
            destinationTetIndices,
            1);
    assert(PHYSIK_IsComponentHandleValid(world, destinationTetMesh) == 1);

    const PhysiK::ComponentHandle mapper =
        PHYSIK_CreateTetMeshMapperComponent(
            world,
            sourceTetMesh,
            destinationTetMesh);
    assert(PHYSIK_IsComponentHandleValid(world, mapper) == 1);

    assert(PHYSIK_SetTetMeshLocalCurrentPosition(
        world,
        sourceTetMesh,
        0,
        0.0f,
        0.0f,
        1.0f) == 1);
    assert(PHYSIK_SetTetMeshLocalCurrentPosition(
        world,
        sourceTetMesh,
        1,
        2.0f,
        0.0f,
        1.0f) == 1);
    assert(PHYSIK_SetTetMeshLocalCurrentPosition(
        world,
        sourceTetMesh,
        2,
        0.0f,
        2.0f,
        1.0f) == 1);
    assert(PHYSIK_SetTetMeshLocalCurrentPosition(
        world,
        sourceTetMesh,
        3,
        0.0f,
        0.0f,
        3.0f) == 1);
    assert(PHYSIK_SetTetMeshLocalCurrentPosition(
        world,
        sourceTetMesh,
        -1,
        0.0f,
        0.0f,
        0.0f) == 0);
    assert(PHYSIK_SetTetMeshLocalCurrentPosition(
        nullptr,
        sourceTetMesh,
        0,
        0.0f,
        0.0f,
        0.0f) == 0);
    assert(PHYSIK_SetTetMeshLocalCurrentPosition(
        world,
        mapper,
        0,
        0.0f,
        0.0f,
        0.0f) == 0);

    PHYSIK_Step(world, 0.0f);

    const Point mapped0 =
        GetTetMeshLocalCurrentPosition(world, destinationTetMesh, 0);
    const Point mapped1 =
        GetTetMeshLocalCurrentPosition(world, destinationTetMesh, 1);
    const Point unmapped =
        GetTetMeshLocalCurrentPosition(world, destinationTetMesh, 4);
    assert(NearlyEqual(mapped0.x, 0.5f));
    assert(NearlyEqual(mapped0.y, 0.5f));
    assert(NearlyEqual(mapped0.z, 1.5f));
    assert(NearlyEqual(mapped1.x, 1.0f));
    assert(NearlyEqual(mapped1.y, 0.5f));
    assert(NearlyEqual(mapped1.z, 1.5f));
    assert(NearlyEqual(unmapped.x, 2.0f));
    assert(NearlyEqual(unmapped.y, 2.0f));
    assert(NearlyEqual(unmapped.z, 2.0f));
    assert(NearlyEqual(
        GetTetMeshLocalCurrentPosition(world, sourceTetMesh, 0).z,
        1.0f));
    Point scratch;
    assert(PHYSIK_GetTetMeshLocalCurrentPosition(
        nullptr,
        destinationTetMesh,
        0,
        &scratch.x,
        &scratch.y,
        &scratch.z) == 0);
    assert(PHYSIK_GetTetMeshLocalCurrentPosition(
        world,
        mapper,
        0,
        &scratch.x,
        &scratch.y,
        &scratch.z) == 0);
    assert(PHYSIK_GetTetMeshLocalCurrentPosition(
        world,
        destinationTetMesh,
        -1,
        &scratch.x,
        &scratch.y,
        &scratch.z) == 0);
    assert(PHYSIK_GetTetMeshLocalCurrentPosition(
        world,
        destinationTetMesh,
        0,
        nullptr,
        &scratch.y,
        &scratch.z) == 0);

    PHYSIK_DestroyWorld(world);
}

void TetMeshMapperComponentRefreshesActiveStatesAfterSourceTopologyChanges()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    const PhysiK::Vec3 sourcePositions[] = {
        PhysiK::Vec3{0.0f, 0.0f, 0.0f},
        PhysiK::Vec3{1.0f, 0.0f, 0.0f},
        PhysiK::Vec3{0.0f, 1.0f, 0.0f},
        PhysiK::Vec3{0.0f, 0.0f, 1.0f},
        PhysiK::Vec3{0.0f, 0.0f, 2.0f}};
    const int sourceTetIndices[] = {
        0, 1, 2, 3,
        0, 1, 2, 4};
    const PhysiK::ComponentHandle sourceTetMesh =
        PHYSIK_CreateTetMeshComponent(
            world,
            sourcePositions,
            5,
            sourceTetIndices,
            2);
    assert(PHYSIK_IsComponentHandleValid(world, sourceTetMesh) == 1);

    const PhysiK::Vec3 destinationPositions[] = {
        PhysiK::Vec3{0.25f, 0.25f, 0.25f},
        PhysiK::Vec3{0.30f, 0.25f, 0.25f},
        PhysiK::Vec3{0.25f, 0.30f, 0.25f},
        PhysiK::Vec3{0.25f, 0.25f, 0.30f}};
    const int destinationTetIndices[] = {0, 1, 2, 3};
    const PhysiK::ComponentHandle destinationTetMesh =
        PHYSIK_CreateTetMeshComponent(
            world,
            destinationPositions,
            4,
            destinationTetIndices,
            1);
    assert(PHYSIK_IsComponentHandleValid(world, destinationTetMesh) == 1);

    const PhysiK::ComponentHandle mapper =
        PHYSIK_CreateTetMeshMapperComponent(
            world,
            sourceTetMesh,
            destinationTetMesh);
    assert(PHYSIK_IsComponentHandleValid(world, mapper) == 1);

    assert(PHYSIK_SetTetMeshLocalCurrentPosition(
        world,
        sourceTetMesh,
        4,
        0.0f,
        0.0f,
        10.0f) == 1);

    PHYSIK_Step(world, 0.0f);
    const Point initiallyMapped =
        GetTetMeshLocalCurrentPosition(world, destinationTetMesh, 0);
    assert(NearlyEqual(initiallyMapped.x, 0.25f));
    assert(NearlyEqual(initiallyMapped.y, 0.25f));
    assert(NearlyEqual(initiallyMapped.z, 0.25f));

    PHYSIK_DeactivateTet(world, sourceTetMesh, 0);
    PHYSIK_Step(world, 0.0f);

    const Point notRemapped =
        GetTetMeshLocalCurrentPosition(world, destinationTetMesh, 0);
    assert(NearlyEqual(notRemapped.x, 0.25f));
    assert(NearlyEqual(notRemapped.y, 0.25f));
    assert(NearlyEqual(notRemapped.z, 0.25f));
    assert(PHYSIK_IsTetActive(world, destinationTetMesh, 0) == 0);

    PHYSIK_DestroyWorld(world);
}

void TetMeshMapperComponentDeactivatesDestinationTetWithInactiveMappedVertex()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    const PhysiK::Vec3 sourcePositions[] = {
        PhysiK::Vec3{0.0f, 0.0f, 0.0f},
        PhysiK::Vec3{0.2f, 0.0f, 0.0f},
        PhysiK::Vec3{0.0f, 0.2f, 0.0f},
        PhysiK::Vec3{0.0f, 0.0f, 0.2f},
        PhysiK::Vec3{3.0f, 0.0f, 0.0f},
        PhysiK::Vec3{0.0f, 3.0f, 0.0f},
        PhysiK::Vec3{0.0f, 0.0f, 3.0f}};
    const int sourceTetIndices[] = {
        0, 1, 2, 3,
        0, 4, 5, 6};
    const PhysiK::ComponentHandle sourceTetMesh =
        PHYSIK_CreateTetMeshComponent(
            world,
            sourcePositions,
            7,
            sourceTetIndices,
            2);
    assert(PHYSIK_IsComponentHandleValid(world, sourceTetMesh) == 1);

    const PhysiK::Vec3 destinationPositions[] = {
        PhysiK::Vec3{0.05f, 0.05f, 0.05f},
        PhysiK::Vec3{1.0f, 0.1f, 0.1f},
        PhysiK::Vec3{0.1f, 1.0f, 0.1f},
        PhysiK::Vec3{0.1f, 0.1f, 1.0f}};
    const int destinationTetIndices[] = {0, 1, 2, 3};
    const PhysiK::ComponentHandle destinationTetMesh =
        PHYSIK_CreateTetMeshComponent(
            world,
            destinationPositions,
            4,
            destinationTetIndices,
            1);
    assert(PHYSIK_IsComponentHandleValid(world, destinationTetMesh) == 1);

    const PhysiK::ComponentHandle mapper =
        PHYSIK_CreateTetMeshMapperComponent(
            world,
            sourceTetMesh,
            destinationTetMesh);
    assert(PHYSIK_IsComponentHandleValid(world, mapper) == 1);

    PHYSIK_Step(world, 0.0f);
    assert(PHYSIK_IsTetActive(world, destinationTetMesh, 0) == 1);

    PHYSIK_DeactivateTet(world, sourceTetMesh, 0);
    PHYSIK_Step(world, 0.0f);

    assert(PHYSIK_IsTetActive(world, destinationTetMesh, 0) == 0);

    PHYSIK_DestroyWorld(world);
}

void TetMeshMapperComponentPropagatesSourceCutsToDestinationTets()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    const PhysiK::Vec3 sourcePositions[] = {
        PhysiK::Vec3{0.0f, 0.0f, 0.0f},
        PhysiK::Vec3{1.0f, 0.0f, 0.0f},
        PhysiK::Vec3{0.0f, 1.0f, 0.0f},
        PhysiK::Vec3{0.0f, 0.0f, 1.0f},
        PhysiK::Vec3{2.0f, 0.0f, 0.0f},
        PhysiK::Vec3{3.0f, 0.0f, 0.0f},
        PhysiK::Vec3{2.0f, 1.0f, 0.0f},
        PhysiK::Vec3{2.0f, 0.0f, 1.0f}};
    const int sourceTetIndices[] = {
        0, 1, 2, 3,
        4, 5, 6, 7};
    const PhysiK::ComponentHandle sourceTetMesh =
        PHYSIK_CreateTetMeshComponent(
            world,
            sourcePositions,
            8,
            sourceTetIndices,
            2);
    assert(PHYSIK_IsComponentHandleValid(world, sourceTetMesh) == 1);

    const PhysiK::Vec3 destinationPositions[] = {
        PhysiK::Vec3{0.10f, 0.10f, 0.10f},
        PhysiK::Vec3{0.40f, 0.10f, 0.10f},
        PhysiK::Vec3{0.10f, 0.40f, 0.10f},
        PhysiK::Vec3{0.10f, 0.10f, 0.40f},
        PhysiK::Vec3{2.10f, 0.10f, 0.10f},
        PhysiK::Vec3{2.40f, 0.10f, 0.10f},
        PhysiK::Vec3{2.10f, 0.40f, 0.10f},
        PhysiK::Vec3{2.10f, 0.10f, 0.40f}};
    const int destinationTetIndices[] = {
        0, 1, 2, 3,
        4, 5, 6, 7};
    const PhysiK::ComponentHandle destinationTetMesh =
        PHYSIK_CreateTetMeshComponent(
            world,
            destinationPositions,
            8,
            destinationTetIndices,
            2);
    assert(PHYSIK_IsComponentHandleValid(world, destinationTetMesh) == 1);

    const PhysiK::ComponentHandle mapper =
        PHYSIK_CreateTetMeshMapperComponent(
            world,
            sourceTetMesh,
            destinationTetMesh);
    assert(PHYSIK_IsComponentHandleValid(world, mapper) == 1);

    PHYSIK_Step(world, 0.0f);
    assert(PHYSIK_IsTetActive(world, destinationTetMesh, 0) == 1);
    assert(PHYSIK_IsTetActive(world, destinationTetMesh, 1) == 1);

    PHYSIK_DeactivateTet(world, sourceTetMesh, 0);
    PHYSIK_Step(world, 0.0f);

    assert(PHYSIK_IsTetActive(world, destinationTetMesh, 0) == 0);
    assert(PHYSIK_IsTetActive(world, destinationTetMesh, 1) == 1);

    PHYSIK_DestroyWorld(world);
}

void TetMeshMapperComponentCanBeCreatedThroughNativeApi()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    int sourceNodes[4] = {};
    const PhysiK::ComponentHandle sourceTetMesh =
        CreateSingleTetMesh(world, sourceNodes);
    const int destinationNodes[] = {
        AddNode(world, 0.25f, 0.25f, 0.25f),
        AddNode(world, 0.50f, 0.25f, 0.25f),
        AddNode(world, 0.25f, 0.50f, 0.25f),
        AddNode(world, 0.25f, 0.25f, 0.50f)};
    const int destinationTetIndices[] = {
        destinationNodes[0],
        destinationNodes[1],
        destinationNodes[2],
        destinationNodes[3]};
    PhysikMaterialDesc material = MakeMaterialDesc(1.0f, 0.0f);
    const PhysiK::ComponentHandle destinationTetMesh =
        PHYSIK_CreateTetMeshPhysicsComponent(
            world,
            destinationNodes,
            4,
            destinationTetIndices,
            1,
            &material,
            0);
    assert(PHYSIK_IsComponentHandleValid(world, destinationTetMesh) == 1);

    const PhysiK::ComponentHandle mapper =
        PHYSIK_CreateTetMeshMapperComponent(
            world,
            sourceTetMesh,
            destinationTetMesh);
    assert(PHYSIK_IsComponentHandleValid(world, mapper) == 1);
    PHYSIK_Step(world, 0.0f);
    assert(PHYSIK_CreateTetMeshMapperComponent(
        nullptr,
        sourceTetMesh,
        destinationTetMesh).IsValid() == false);
    assert(PHYSIK_CreateTetMeshMapperComponent(
        world,
        PhysiK::ComponentHandle{},
        destinationTetMesh).IsValid() == false);

    PHYSIK_DestroyWorld(world);
}

void FemModelLinearRouteUsesExistingAssembly()
{
    std::vector<PhysiK::Node> nodes = CreateUnitTetNodes();
    PhysiK::Tet tet = CreateUnitTet();
    PhysiK::FEMModel::InitializeTetRestData(tet, RestPositionsFromNodes(nodes));
    nodes[3].position = PhysiK::Vec3{0.0f, 0.0f, 1.1f};

    PhysiK::SolverData solverData;
    const bool implemented = ComputeFemForcesIntoSolverData(
        PhysiK::FemModel::Linear,
        {tet},
        nodes,
        solverData);

    assert(implemented);
    assert(!solverData.GetNodeForces().empty());
    assert(!solverData.GetStiffnessBlocks().empty());
    assert(SumForcesForNode(solverData, 3).z < 0.0f);
}

void FemModelCorotationalRouteUsesCorotationalAssembly()
{
    std::vector<PhysiK::Node> nodes = CreateUnitTetNodes();
    PhysiK::Tet tet = CreateUnitTet();
    PhysiK::FEMModel::InitializeTetRestData(tet, RestPositionsFromNodes(nodes));
    nodes[3].position = PhysiK::Vec3{0.0f, 0.0f, 1.1f};

    PhysiK::SolverData solverData;
    const bool implemented = ComputeFemForcesIntoSolverData(
        PhysiK::FemModel::Corotational,
        {tet},
        nodes,
        solverData);

    assert(implemented);
    assert(PhysiK::FEMModel::IsFemModelImplemented(PhysiK::FemModel::Corotational));
    assert(!solverData.GetNodeForces().empty());
    assert(!solverData.GetStiffnessBlocks().empty());
    assert(SumForcesForNode(solverData, 3).z < 0.0f);
}

void FemModelNeoHookeanRouteIsExplicitlyNotImplemented()
{
    std::vector<PhysiK::Node> nodes = CreateUnitTetNodes();
    PhysiK::Tet tet = CreateUnitTet();
    PhysiK::FEMModel::InitializeTetRestData(tet, RestPositionsFromNodes(nodes));
    nodes[3].position = PhysiK::Vec3{0.0f, 0.0f, 1.1f};

    PhysiK::SolverData solverData;
    const bool implemented = ComputeFemForcesIntoSolverData(
        PhysiK::FemModel::NeoHookean,
        {tet},
        nodes,
        solverData);

    assert(!implemented);
    assert(!PhysiK::FEMModel::IsFemModelImplemented(PhysiK::FemModel::NeoHookean));
    assert(solverData.GetNodeForces().empty());
    assert(solverData.GetStiffnessBlocks().empty());
}

void CorotationalFemHasNearZeroForceForRigidRotation()
{
    std::vector<PhysiK::Node> nodes = CreateUnitTetNodes();
    PhysiK::Tet tet = CreateUnitTet(1000.0f);
    PhysiK::FEMModel::InitializeTetRestData(tet, RestPositionsFromNodes(nodes));

    for (PhysiK::Node& node : nodes)
    {
        node.position = RotateZ90(node.position);
    }

    PhysiK::SolverData solverData;
    const bool implemented = ComputeFemForcesIntoSolverData(
        PhysiK::FemModel::Corotational,
        {tet},
        nodes,
        solverData);

    assert(implemented);
    assert(SumForceLengthSquared(solverData) < 0.000001f);
}

void LinearFemProducesForceForRigidRotation()
{
    std::vector<PhysiK::Node> nodes = CreateUnitTetNodes();
    PhysiK::Tet tet = CreateUnitTet(1000.0f);
    PhysiK::FEMModel::InitializeTetRestData(tet, RestPositionsFromNodes(nodes));

    for (PhysiK::Node& node : nodes)
    {
        node.position = RotateZ90(node.position);
    }

    PhysiK::SolverData solverData;
    const bool implemented = ComputeFemForcesIntoSolverData(
        PhysiK::FemModel::Linear,
        {tet},
        nodes,
        solverData);

    assert(implemented);
    assert(SumForceLengthSquared(solverData) > 0.000001f);
}

void CorotationalFemProducesRestoringForceForSmallDeformation()
{
    std::vector<PhysiK::Node> nodes = CreateUnitTetNodes();
    PhysiK::Tet tet = CreateUnitTet(1000.0f);
    PhysiK::FEMModel::InitializeTetRestData(tet, RestPositionsFromNodes(nodes));
    nodes[3].position.z += 0.1f;

    PhysiK::SolverData solverData;
    const bool implemented = ComputeFemForcesIntoSolverData(
        PhysiK::FemModel::Corotational,
        {tet},
        nodes,
        solverData);

    assert(implemented);
    const PhysiK::Vec3 node3Force = SumForcesForNode(solverData, 3);
    assert(LengthSquared(node3Force) > 0.000001f);
    assert(node3Force.z < 0.0f);
}

void CorotationalAssemblySmokeTestStaysFinite()
{
    std::vector<PhysiK::Node> nodes = CreateUnitTetNodes();
    PhysiK::Tet tet = CreateUnitTet(250.0f);
    PhysiK::FEMModel::InitializeTetRestData(tet, RestPositionsFromNodes(nodes));

    for (PhysiK::Node& node : nodes)
    {
        node.position = RotateZ90(node.position);
    }
    nodes[3].position.z += 0.05f;

    PhysiK::SolverData solverData;
    const bool implemented = ComputeFemForcesIntoSolverData(
        PhysiK::FemModel::Corotational,
        {tet},
        nodes,
        solverData);

    assert(implemented);
    assert(!solverData.GetNodeForces().empty());
    assert(!solverData.GetStiffnessBlocks().empty());
    for (const PhysiK::SolverData::NodeForce& force : solverData.GetNodeForces())
    {
        assert(IsFinite(force.force));
    }
    for (const PhysiK::SolverData::StiffnessBlock& block : solverData.GetStiffnessBlocks())
    {
        assert(IsFinite(block.block));
    }
}

void LinearTetAssemblyProducesForcesAndSymmetricStiffness()
{
    std::vector<PhysiK::Node> nodes = CreateUnitTetNodes();
    PhysiK::Tet tet = CreateUnitTet();
    PhysiK::FEMModel::InitializeTetRestData(tet, RestPositionsFromNodes(nodes));

    std::vector<PhysiK::Tet> tets = {tet};
    PhysiK::SolverData solverData;
    ComputeLinearFemForcesIntoSolverData(tets, nodes, solverData);

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
    ComputeLinearFemForcesIntoSolverData(tets, nodes, solverData);

    const PhysiK::Vec3 node3Force = SumForcesForNode(solverData, 3);
    assert(LengthSquared(node3Force) > 0.000001f);
    assert(node3Force.z < 0.0f);
}

void UnitTetShapeFunctionGradientsMatchExpectedConvention()
{
    std::vector<PhysiK::Node> nodes = CreateUnitTetNodes();
    PhysiK::Tet tet = CreateUnitTet();
    PhysiK::FEMModel::InitializeTetRestData(tet, RestPositionsFromNodes(nodes));

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
    PhysiK::FEMModel::InitializeTetRestData(tet, RestPositionsFromNodes(nodes));

    assert(tet.restVolume == 0.0f);
    for (const PhysiK::Vec3& gradient : tet.shapeFunctionGradients)
    {
        assert(NearlyEqual(gradient, PhysiK::Vec3{}));
    }

    std::vector<PhysiK::Tet> tets = {tet};
    PhysiK::SolverData solverData;
    ComputeLinearFemForcesIntoSolverData(tets, nodes, solverData);

    assert(solverData.GetNodeForces().empty());
    assert(solverData.GetStiffnessBlocks().empty());
}

void LinearTetMaterialSanitizationAvoidsInvalidForces()
{
    std::vector<PhysiK::Node> nodes = CreateUnitTetNodes();
    PhysiK::Tet noElasticity = CreateUnitTet(-10.0f, 0.25f);
    PhysiK::FEMModel::InitializeTetRestData(noElasticity, RestPositionsFromNodes(nodes));
    nodes[1].position.x += 0.1f;

    PhysiK::SolverData solverData;
    ComputeLinearFemForcesIntoSolverData({noElasticity}, nodes, solverData);

    for (const PhysiK::SolverData::NodeForce& force : solverData.GetNodeForces())
    {
        assert(IsFinite(force.force));
        assert(LengthSquared(force.force) < 0.000001f);
    }

    nodes = CreateUnitTetNodes();
    PhysiK::Tet clampedPoisson = CreateUnitTet(100.0f, 0.99f);
    PhysiK::FEMModel::InitializeTetRestData(clampedPoisson, RestPositionsFromNodes(nodes));
    nodes[1].position.x += 0.1f;
    solverData.Clear();
    ComputeLinearFemForcesIntoSolverData({clampedPoisson}, nodes, solverData);

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
    PhysiK::FEMModel::InitializeTetRestData(tet, RestPositionsFromNodes(nodes));
    nodes[1].position.x += 0.1f;

    PhysiK::SolverData solverData;
    ComputeLinearFemForcesIntoSolverData({tet}, nodes, solverData);

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
    SetNodeFixedControlsDynamicState();
    SphereContactCreatesTransientConnectionAndMovesTet();
    CollisionSphereConnectionSettingsHaveDefaults();
    CollisionSphereConnectionSettingsCanBeUpdatedAndClamped();
    CollisionSphereConnectionSettingsAffectGeneratedConnections();
    CollisionSphereConnectionSettingsCAPIHandlesInvalidInputs();
    MultipleForceSourcesCoexist();
    ExternalLogicHookRunsOnceBeforeSubsteps();
    KinematicUpdateRunsAfterExternalLogicBeforePhysicsSubsteps();
    FEMElasticityMovesDistortedTetTowardRestShape();
    FEMUsesWorldNodeMappingWhenLocalTetIndicesDiffer();
    TetMeshPhysicsUsesCreationRestDataForRuntimeForces();
    ImplicitEulerFEMTetMovesDistortedNodeTowardRestShape();
    ImplicitEulerAllDynamicTetResistsRestShapeVelocity();
    ImplicitEulerPointAnchoredTetRecoversFreeNodes();
    ImplicitEulerFEMRecoveryBeatsGravityOnlyMotion();
    FEMLumpedMassUsesDensityAndRestVolume();
    TetMeshMaterialDescCrossesNativeBoundary();
    TetMeshMaterialCanBeUpdatedThroughNativeDescriptor();
    FEMLumpedMassDensityReducesAccelerationFromSameForce();
    FEMLumpedMassPreservesFixedNodes();
    FEMGravityAccelerationIsIndependentOfDensity();
    ImplicitEulerUsesStiffnessBlocks();
    ImplicitPointConnectionStableWithHighStiffness();
    PointConnectionBarycentricAssemblyDistributesForcesAndStiffness();
    ImplicitAnchoredTetPointConnectionsRemainStableUnderGravity();
    TetMeshComponentOwnsTetsAndWorldStepUsesComponentSystem();
    CollisionSphereOverlapQueryOverlapsOneTetNode();
    CollisionSphereOverlapUsesWorldNodeMapping();
    CollisionSphereOverlapQueryReturnsZeroWhenNoNodesOverlap();
    CollisionSphereOverlapQueryIgnoresInactiveTets();
    CollisionSphereOverlapQueryIgnoresDestroyedTetMesh();
    CollisionSphereOverlapQueryIgnoresZeroRadiusSphere();
    CollisionSphereOverlapQueryReportsMaskAndCount();
    CollisionSphereOverlapQueryCountAndFillAreConsistent();
    CollisionSphereOverlapQueryInvalidInputsReturnZero();
    SparseBlockMatrixStoresAndMultipliesBlocks();
    SparseBlockMatrixAddBlockAccumulatesContributions();
    SparseBlockMatrixSingleTetPatternContainsAllCouplings();
    SparseBlockMatrixAdjacentTetsReuseSharedBlocks();
    TetMeshComponentCachesFemSparsePattern();
    TetMeshComponentMapsLocalFemPatternToGlobalSolverNodes();
    ConjugateGradientSolvesDiagonalSparseSystem();
    ConjugateGradientSolvesCoupledSparseSystem();
    CurrentLinearSolverSolvesKnownSparseSystem();
    SolverDataFailedImplicitSolveLeavesNoDeltaVelocity();
    PerformanceLoggingWritesCsvForImplicitStep();
    ImplicitEulerLinearTetUsesSparseCgPath();
    ImplicitEulerCorotationalTetUsesSparseCgPath();
    MultiTetImplicitEulerSparseCgSmokeTest();
    TetActiveStateDefaultsAndNoOpsAreSafe();
    TetActiveStateIsExposedThroughNativeApi();
    InactiveTetsAreSkippedByFemForceAndStiffnessAssembly();
    DeactivatedTetsAreSkippedByLumpedMassAssembly();
    DeactivatingTetDoesNotDirtySparsePattern();
    SmallTetMeshSimulatesAfterTetDeactivation();
    TetMeshComponentStoresGeometryWithoutWorldNodes();
    TetMeshComponentDefaultFemModelIsLinear();
    TetMeshComponentStoresSelectedFemModel();
    EventSystemDeliversSubscribedEventsOnlyOnce();
    TopologyMeshComponentDeclaresEventsAndClearsDirtyFlag();
    TopologyMeshComponentBuildsActiveTetIslands();
    SurfaceExtractionComponentExtractsActiveBoundaryFaces();
    SurfaceExtractionComponentWindsBoundaryFacesOutward();
    TetMeshGeneratorWeldsNodesAndDropsDegenerateTets();
    GeneratedTetMeshApiStoresCleanGeometry();
    TetMeshComponentCanBeCreatedFromGeneratedTetMeshApi();
    TetMeshPhysicsComponentCanBeCreatedFromGeneratedTetMeshApi();
    GeneratedTetMeshApiWeldsDuplicateNodes();
    TetMeshPhysicsCanConsumeGeneratedTetMeshDirectly();
    TetMeshCreationWeldsDuplicateSharedFaceNodes();
    SurfaceExtractionComponentCanBeCreatedThroughNativeApi();
    SurfaceVisualComponentBuildsRenderReadySurface();
    VisualMeshComponentDeclaresTopologyListenerAndClearsDirtyFlag();
    VisualMeshComponentCanBeCreatedThroughNativeApi();
    VisualMeshComponentStoresVisualMeshData();
    VisualMeshComponentBuildsBruteForceEmbedding();
    VisualMeshComponentUpdatesDeformedVerticesFromHostTet();
    VisualMeshComponentCAPIExportsMeshBuffers();
    TetMeshMapperComponentEmbedsDestinationAndFollowsSource();
    TetMeshMapperComponentRefreshesActiveStatesAfterSourceTopologyChanges();
    TetMeshMapperComponentDeactivatesDestinationTetWithInactiveMappedVertex();
    TetMeshMapperComponentPropagatesSourceCutsToDestinationTets();
    TetMeshMapperComponentCanBeCreatedThroughNativeApi();
    FemModelLinearRouteUsesExistingAssembly();
    FemModelCorotationalRouteUsesCorotationalAssembly();
    FemModelNeoHookeanRouteIsExplicitlyNotImplemented();
    CorotationalFemHasNearZeroForceForRigidRotation();
    LinearFemProducesForceForRigidRotation();
    CorotationalFemProducesRestoringForceForSmallDeformation();
    CorotationalAssemblySmokeTestStaysFinite();
    UnitTetShapeFunctionGradientsMatchExpectedConvention();
    DegenerateTetIsSkippedWithoutInvalidAssembly();
    LinearTetMaterialSanitizationAvoidsInvalidForces();
    LinearTetDisplacedNodeReceivesRestoringForce();
    LinearTetAssemblyProducesForcesAndSymmetricStiffness();
    DestroyComponentInvalidatesHandle();
    return 0;
}
