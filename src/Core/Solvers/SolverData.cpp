#include "PhysiK/Core/Solvers/SolverData.h"

#include <algorithm>
#include <cassert>
#include <cmath>

#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
#include <chrono>
#define PHYSIK_COLLECT_CG_TIMING 1
#endif

namespace PhysiK
{
    namespace
    {
#if defined(PHYSIK_COLLECT_CG_TIMING)
        using Clock = std::chrono::steady_clock;

        double ElapsedMilliseconds(Clock::time_point start)
        {
            return std::chrono::duration<double, std::milli>(
                Clock::now() - start).count();
        }
#endif

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
        ClearTransientState();
        matrix.Clear();
        implicitPatternDirty = true;
        cachedDynamicBlockCount = 0;
        cachedNodeToDynamicBlock.clear();
        cachedBlockCoordinates.clear();
    }

    void SolverData::ClearTransientState()
    {
        nodeForces.clear();
        nodeMasses.clear();
        stiffnessBlocks.clear();
        assembledMasses.clear();
        assembledForces.clear();
        nodeToDynamicBlock.clear();
        dynamicBlockCount = 0;
        rhs.clear();
        deltaVelocity.clear();
        inversePreconditioner.clear();
        cgResidual.clear();
        cgDirection.clear();
        cgTemp.clear();
        lastLinearSolveResult = LinearSolveResult{};
    }

    void SolverData::MarkImplicitPatternDirty()
    {
        implicitPatternDirty = true;
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

    bool SolverData::PrecomputeImplicitSolve(
        const std::vector<Node>& nodes,
        const std::vector<Vec3>& nodeVelocities,
        float dt)
    {
        // World::BuildSolverData assembles masses before validation and gravity.
        // Keep this fallback so direct SolverData users and tests remain safe.
        if (assembledMasses.size() != nodes.size())
        {
            AssembleMasses(static_cast<int>(nodes.size()));
        }

        assembledForces.assign(nodes.size(), Vec3{});
        nodeToDynamicBlock.clear();
        dynamicBlockCount = 0;
        rhs.clear();
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

        return AssembleImplicitMatrixAndRhs(nodes, nodeVelocities, dt);
    }

    bool SolverData::SolveImplicitLinearSystem(const ConjugateGradientSettings& cgSettings)
    {
        if (dynamicBlockCount <= 0)
        {
            return false;
        }

#if defined(PHYSIK_COLLECT_CG_TIMING)
        Clock::time_point timerStart = Clock::now();
#endif
        if (!BuildInversePreconditioner(cgSettings.useJacobiPreconditioner))
        {
            lastLinearSolveResult = LinearSolveResult{};
#if defined(PHYSIK_COLLECT_CG_TIMING)
            lastLinearSolveResult.preconditionerBuildMs =
                ElapsedMilliseconds(timerStart);
#endif
            return false;
        }
        lastLinearSolveResult = LinearSolveResult{};
#if defined(PHYSIK_COLLECT_CG_TIMING)
        lastLinearSolveResult.preconditionerBuildMs =
            ElapsedMilliseconds(timerStart);

        timerStart = Clock::now();
#endif
        const ConjugateGradientResult result =
            SolvePreconditionedConjugateGradient(
                deltaVelocity,
                matrix,
                rhs,
                cgSettings.maxIterations,
                cgSettings.tolerance,
                inversePreconditioner,
                cgResidual,
                cgDirection,
                cgTemp);
#if defined(PHYSIK_COLLECT_CG_TIMING)
        lastLinearSolveResult.cgTotalMs = ElapsedMilliseconds(timerStart);
#endif

        lastLinearSolveResult.iterations = result.iterations;
        lastLinearSolveResult.residualNorm = result.residualNorm;
        lastLinearSolveResult.converged = result.converged;
#if defined(PHYSIK_COLLECT_CG_TIMING)
        lastLinearSolveResult.cgMultiplyMs = result.cgMultiplyMs;
        lastLinearSolveResult.cgApplyPreconditionerMs =
            result.cgApplyPreconditionerMs;
        lastLinearSolveResult.cgDotVectorOpsMs = result.cgDotVectorOpsMs;
#endif

        if (deltaVelocity.size() != static_cast<std::size_t>(dynamicBlockCount))
        {
            return false;
        }

        for (const Vec3& value : deltaVelocity)
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

    bool SolverData::BuildInversePreconditioner(bool useJacobiPreconditioner)
    {
        inversePreconditioner.assign(
            static_cast<std::size_t>(std::max(0, dynamicBlockCount)),
            Mat3::Identity());

        if (!useJacobiPreconditioner)
        {
            return true;
        }

        constexpr float DeterminantTolerance = 1.0e-8f;
        for (int block = 0; block < dynamicBlockCount; ++block)
        {
            const int blockIndex = matrix.FindBlockIndex(block, block);
            if (blockIndex < 0)
            {
                assert(false && "implicit matrix is missing a diagonal block");
                return false;
            }

            const Mat3& diagonalBlock =
                matrix.values[static_cast<std::size_t>(blockIndex)];

            if (!IsFinite(diagonalBlock))
            {
                return false;
            }

            const float determinant = Determinant(diagonalBlock);
            if (!IsFinite(determinant) ||
                std::abs(determinant) <= DeterminantTolerance)
            {
                return false;
            }

            const Mat3 inverseBlock = Inverse(diagonalBlock);
            if (!IsFinite(inverseBlock))
            {
                return false;
            }

            inversePreconditioner[static_cast<std::size_t>(block)] = inverseBlock;
        }

        return true;
    }

    bool SolverData::AssembleImplicitMatrixAndRhs(
        const std::vector<Node>& nodes,
        const std::vector<Vec3>& nodeVelocities,
        float dt)
    {
        rhs.assign(static_cast<std::size_t>(dynamicBlockCount), Vec3{});

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

            if (block.nodeB >= static_cast<int>(nodeVelocities.size()))
            {
                return false;
            }

            const Vec3& columnVelocity = nodeVelocities[static_cast<std::size_t>(block.nodeB)];
            if (!IsFinite(columnVelocity))
            {
                return false;
            }

            const Vec3 stiffnessVelocity = block.block * columnVelocity;
            if (!IsFinite(stiffnessVelocity))
            {
                return false;
            }

            rhs[static_cast<std::size_t>(rowBlock)] -=
                stiffnessVelocity * stiffnessScale;

            const int columnBlock = GetDynamicBlockForNode(block.nodeB);
            if (columnBlock < 0)
            {
                continue;
            }

            blockCoordinates.push_back({rowBlock, columnBlock});
        }

        const bool rebuildPattern = implicitPatternDirty ||
            cachedDynamicBlockCount != dynamicBlockCount ||
            cachedNodeToDynamicBlock != nodeToDynamicBlock ||
            cachedBlockCoordinates != blockCoordinates;

        if (rebuildPattern)
        {
            matrix.BuildPattern(dynamicBlockCount, blockCoordinates);
            cachedDynamicBlockCount = dynamicBlockCount;
            cachedNodeToDynamicBlock = nodeToDynamicBlock;
            cachedBlockCoordinates = blockCoordinates;
            implicitPatternDirty = false;
            ++implicitPatternRebuildCount;
        }
        else
        {
            matrix.ClearValues();
            ++implicitPatternReuseCount;
        }

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
            rhs[static_cast<std::size_t>(dynamicBlock)] += force * dt;
        }

        return true;
    }
}
