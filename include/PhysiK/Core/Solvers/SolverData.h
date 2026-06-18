#pragma once

#include <utility>
#include <vector>

#include "PhysiK/API/PhysiKAPI.h"
#include "PhysiK/Core/Solvers/Linear/LinearSolver.h"
#include "PhysiK/Math/Mat3.h"
#include "PhysiK/Math/SparseBlockMatrix.h"
#include "PhysiK/Math/Vec3.h"
#include "PhysiK/PhysicsData/Node.h"

namespace PhysiK
{
    class PHYSIK_API SolverData
    {
    public:
        struct NodeForce
        {
            int node = -1;
            Vec3 force;
        };

        struct NodeMass
        {
            int node = -1;
            float mass = 0.0f;
        };

        struct StiffnessBlock
        {
            int nodeA = -1;
            int nodeB = -1;
            Mat3 block;
        };

        void Clear();
        void ClearTransientState();
        void MarkImplicitPatternDirty();

        void AddNodeForce(int node, const Vec3& force)
        {
            nodeForces.push_back(NodeForce{node, force});
        }

        void AddNodeMass(int node, float mass)
        {
            nodeMasses.push_back(NodeMass{node, mass});
        }

        void AddStiffnessBlock(int nodeA, int nodeB, const Mat3& block)
        {
            stiffnessBlocks.push_back(StiffnessBlock{nodeA, nodeB, block});
        }

        const std::vector<NodeForce>& GetNodeForces() const
        {
            return nodeForces;
        }

        const std::vector<NodeMass>& GetNodeMasses() const
        {
            return nodeMasses;
        }

        const std::vector<StiffnessBlock>& GetStiffnessBlocks() const
        {
            return stiffnessBlocks;
        }

        void AssembleMasses(int nodeCount);
        bool HasNodeMassContribution(int nodeIndex) const;
        float GetAssembledMassForNode(int nodeIndex) const;

        bool PrecomputeImplicitSolve(
            const std::vector<Node>& nodes,
            const std::vector<Vec3>& nodeVelocities,
            float dt);
        bool SolveImplicitLinearSystem(const ConjugateGradientSettings& settings);
        const LinearSolveResult& GetLastLinearSolveResult() const;
        int GetLastCgIterationCount() const;
        float GetLastCgResidualNorm() const;

        int GetDynamicBlockForNode(int nodeIndex) const;
        int GetDynamicBlockCount() const
        {
            return dynamicBlockCount;
        }

        const std::vector<Vec3>& GetDeltaVelocity() const
        {
            return deltaVelocity;
        }

        const std::vector<float>& GetAssembledMasses() const
        {
            return assembledMasses;
        }

        const std::vector<Vec3>& GetAssembledForces() const
        {
            return assembledForces;
        }

#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
        int GetImplicitPatternRebuildCount() const
        {
            return implicitPatternRebuildCount;
        }

        int GetImplicitPatternReuseCount() const
        {
            return implicitPatternReuseCount;
        }
#endif

    private:
        bool BuildDynamicNodeMapping(const std::vector<Node>& nodes);
        void BuildInversePreconditioner(bool useJacobiPreconditioner);
        bool AssembleImplicitMatrixAndRhs(
            const std::vector<Node>& nodes,
            const std::vector<Vec3>& nodeVelocities,
            float dt);

        std::vector<NodeForce> nodeForces;
        std::vector<NodeMass> nodeMasses;
        std::vector<StiffnessBlock> stiffnessBlocks;

        std::vector<float> assembledMasses;
        std::vector<Vec3> assembledForces;
        std::vector<int> nodeToDynamicBlock;
        int dynamicBlockCount = 0;
        std::vector<Vec3> rhs;
        SparseBlockMatrix matrix;
        std::vector<Vec3> deltaVelocity;
        std::vector<Mat3> inversePreconditioner;
        std::vector<Vec3> cgResidual;
        std::vector<Vec3> cgDirection;
        std::vector<Vec3> cgTemp;
        LinearSolveResult lastLinearSolveResult;

        bool implicitPatternDirty = true;
        int cachedDynamicBlockCount = 0;
        std::vector<int> cachedNodeToDynamicBlock;
        std::vector<std::pair<int, int>> cachedBlockCoordinates;
        int implicitPatternRebuildCount = 0;
        int implicitPatternReuseCount = 0;
    };
}
