#include "PhysiK/Core/Solvers/SolverData.h"

#include <algorithm>
#include <cmath>

namespace PhysiK
{
    namespace
    {
        bool IsFinite(float value)
        {
            return std::isfinite(value);
        }

        bool IsFinite(const Vec3& value)
        {
            return IsFinite(value.x) && IsFinite(value.y) && IsFinite(value.z);
        }

        bool IsFinite(const Mat3& matrix)
        {
            return IsFinite(matrix.columns[0]) &&
                IsFinite(matrix.columns[1]) &&
                IsFinite(matrix.columns[2]);
        }

        Mat3 ScaledIdentity(float value)
        {
            return Mat3::FromColumns(
                Vec3{value, 0.0f, 0.0f},
                Vec3{0.0f, value, 0.0f},
                Vec3{0.0f, 0.0f, value});
        }

        Mat3 Scale(const Mat3& matrix, float value)
        {
            return Mat3::FromColumns(
                matrix.columns[0] * value,
                matrix.columns[1] * value,
                matrix.columns[2] * value);
        }

    }

    void SolverData::Clear()
    {
        nodeForces.clear();
        nodeMasses.clear();
        stiffnessBlocks.clear();
        assembledMasses.clear();
        assembledForces.clear();
        nodeToDynamicBlock.clear();
        dynamicBlockCount = 0;
        rhs.clear();
        matrix.Clear();
        deltaVelocity.clear();
        lastLinearSolveResult = LinearSolveResult{};
    }

    void SolverData::AssembleMasses(int nodeCount)
    {
        assembledMasses.assign(static_cast<std::size_t>(std::max(0, nodeCount)), 0.0f);

        for (const NodeMass& nodeMass : nodeMasses)
        {
            if (nodeMass.node < 0 || nodeMass.node >= nodeCount ||
                !std::isfinite(nodeMass.mass) || nodeMass.mass <= 0.0f)
            {
                continue;
            }

            assembledMasses[static_cast<std::size_t>(nodeMass.node)] += nodeMass.mass;
        }
    }

    bool SolverData::HasNodeMassContribution(int nodeIndex) const
    {
        for (const NodeMass& nodeMass : nodeMasses)
        {
            if (nodeMass.node == nodeIndex)
            {
                return true;
            }
        }

        return false;
    }

    float SolverData::GetAssembledMassForNode(int nodeIndex) const
    {
        if (nodeIndex < 0 || nodeIndex >= static_cast<int>(assembledMasses.size()))
        {
            return 0.0f;
        }

        return assembledMasses[static_cast<std::size_t>(nodeIndex)];
    }

    bool SolverData::PrecomputeImplicitSolve(const std::vector<Node>& nodes, float dt)
    {
        AssembleMasses(static_cast<int>(nodes.size()));
        assembledForces.assign(nodes.size(), Vec3{});
        nodeToDynamicBlock.clear();
        dynamicBlockCount = 0;
        rhs.clear();
        matrix.Clear();
        deltaVelocity.clear();
        lastLinearSolveResult = LinearSolveResult{};

        for (const NodeForce& nodeForce : nodeForces)
        {
            if (nodeForce.node < 0 || nodeForce.node >= static_cast<int>(nodes.size()) ||
                !IsFinite(nodeForce.force))
            {
                continue;
            }

            assembledForces[static_cast<std::size_t>(nodeForce.node)] += nodeForce.force;
        }

        if (dt <= 0.0f)
        {
            return false;
        }

        if (!BuildDynamicNodeMapping(nodes))
        {
            return false;
        }

        return AssembleImplicitMatrixAndRhs(nodes, dt);
    }

    bool SolverData::SolveImplicitLinearSystem(const ConjugateGradientSettings& cgSettings)
    {
        const std::size_t dimension = static_cast<std::size_t>(dynamicBlockCount * 3);
        if (dimension == 0)
        {
            return false;
        }

        LinearSolveSettings settings;
        settings.maxIterations = cgSettings.maxIterations;
        settings.tolerance = cgSettings.tolerance;
        settings.useJacobiPreconditioner = cgSettings.useJacobiPreconditioner;

        const LinearSolveResult result =
            GetCurrentLinearSolver().Solve(matrix, rhs, deltaVelocity, settings);
        lastLinearSolveResult = result;

        if (deltaVelocity.size() != dimension)
        {
            return false;
        }

        for (float value : deltaVelocity)
        {
            if (!IsFinite(value))
            {
                return false;
            }
        }

        return true;
    }

    const LinearSolveResult& SolverData::GetLastLinearSolveResult() const
    {
        return lastLinearSolveResult;
    }

    int SolverData::GetLastCgIterationCount() const
    {
        return lastLinearSolveResult.iterations;
    }

    float SolverData::GetLastCgResidualNorm() const
    {
        return lastLinearSolveResult.residualNorm;
    }

    int SolverData::GetDynamicBlockForNode(int nodeIndex) const
    {
        if (nodeIndex < 0 || nodeIndex >= static_cast<int>(nodeToDynamicBlock.size()))
        {
            return -1;
        }

        return nodeToDynamicBlock[static_cast<std::size_t>(nodeIndex)];
    }

    bool SolverData::BuildDynamicNodeMapping(const std::vector<Node>& nodes)
    {
        nodeToDynamicBlock.assign(nodes.size(), -1);
        dynamicBlockCount = 0;

        for (int nodeIndex = 0; nodeIndex < static_cast<int>(nodes.size()); ++nodeIndex)
        {
            const float mass = assembledMasses[static_cast<std::size_t>(nodeIndex)];
            const Node& node = nodes[static_cast<std::size_t>(nodeIndex)];
            if (node.active &&
                !node.fixed &&
                std::isfinite(mass) && mass > 0.0f)
            {
                nodeToDynamicBlock[static_cast<std::size_t>(nodeIndex)] = dynamicBlockCount;
                ++dynamicBlockCount;
            }
        }

        return dynamicBlockCount > 0;
    }

    bool SolverData::AssembleImplicitMatrixAndRhs(
        const std::vector<Node>& nodes,
        float dt)
    {
        const std::size_t dimension = static_cast<std::size_t>(dynamicBlockCount * 3);
        rhs.assign(dimension, 0.0f);

        std::vector<std::pair<int, int>> blockCoordinates;
        blockCoordinates.reserve(
            static_cast<std::size_t>(dynamicBlockCount) + stiffnessBlocks.size());

        for (int nodeIndex = 0; nodeIndex < static_cast<int>(nodes.size()); ++nodeIndex)
        {
            const int dynamicBlock = GetDynamicBlockForNode(nodeIndex);
            if (dynamicBlock < 0)
            {
                continue;
            }

            const float mass = assembledMasses[static_cast<std::size_t>(nodeIndex)];
            if (!IsFinite(mass))
            {
                return false;
            }

            blockCoordinates.push_back({dynamicBlock, dynamicBlock});
        }

        const float stiffnessScale = dt * dt;
        for (const StiffnessBlock& block : stiffnessBlocks)
        {
            if (block.nodeA < 0 || block.nodeA >= static_cast<int>(nodes.size()) ||
                block.nodeB < 0 || block.nodeB >= static_cast<int>(nodes.size()))
            {
                continue;
            }

            const int rowBlock = GetDynamicBlockForNode(block.nodeA);
            if (rowBlock < 0)
            {
                continue;
            }

            if (!IsFinite(block.block))
            {
                return false;
            }

            const Vec3& columnVelocity = nodes[static_cast<std::size_t>(block.nodeB)].velocity;
            if (!IsFinite(columnVelocity))
            {
                return false;
            }

            const Vec3 stiffnessVelocity = block.block * columnVelocity;
            if (!IsFinite(stiffnessVelocity))
            {
                return false;
            }

            const int rowBase = rowBlock * 3;
            rhs[static_cast<std::size_t>(rowBase + 0)] -= stiffnessScale * stiffnessVelocity.x;
            rhs[static_cast<std::size_t>(rowBase + 1)] -= stiffnessScale * stiffnessVelocity.y;
            rhs[static_cast<std::size_t>(rowBase + 2)] -= stiffnessScale * stiffnessVelocity.z;

            const int columnBlock = GetDynamicBlockForNode(block.nodeB);
            if (columnBlock < 0)
            {
                continue;
            }

            blockCoordinates.push_back({rowBlock, columnBlock});
        }

        matrix.BuildPattern(dynamicBlockCount, blockCoordinates);

        for (int nodeIndex = 0; nodeIndex < static_cast<int>(nodes.size()); ++nodeIndex)
        {
            const int dynamicBlock = GetDynamicBlockForNode(nodeIndex);
            if (dynamicBlock < 0)
            {
                continue;
            }

            const float mass = assembledMasses[static_cast<std::size_t>(nodeIndex)];
            if (!matrix.AddBlock(dynamicBlock, dynamicBlock, ScaledIdentity(mass)))
            {
                return false;
            }
        }

        for (const StiffnessBlock& block : stiffnessBlocks)
        {
            if (block.nodeA < 0 || block.nodeA >= static_cast<int>(nodes.size()) ||
                block.nodeB < 0 || block.nodeB >= static_cast<int>(nodes.size()) ||
                !IsFinite(block.block))
            {
                continue;
            }

            const int rowBlock = GetDynamicBlockForNode(block.nodeA);
            const int columnBlock = GetDynamicBlockForNode(block.nodeB);
            if (rowBlock < 0 || columnBlock < 0)
            {
                continue;
            }

            if (!matrix.AddBlock(rowBlock, columnBlock, Scale(block.block, stiffnessScale)))
            {
                return false;
            }
        }

        for (int nodeIndex = 0; nodeIndex < static_cast<int>(nodes.size()); ++nodeIndex)
        {
            const int dynamicBlock = GetDynamicBlockForNode(nodeIndex);
            if (dynamicBlock < 0)
            {
                continue;
            }

            const Vec3& force = assembledForces[static_cast<std::size_t>(nodeIndex)];
            const int baseDof = dynamicBlock * 3;
            rhs[static_cast<std::size_t>(baseDof + 0)] += dt * force.x;
            rhs[static_cast<std::size_t>(baseDof + 1)] += dt * force.y;
            rhs[static_cast<std::size_t>(baseDof + 2)] += dt * force.z;
        }

        return true;
    }
}
