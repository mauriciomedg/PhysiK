#pragma once

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

        bool PrecomputeImplicitSolve(const std::vector<Node>& nodes, float dt);
        bool SolveImplicitLinearSystem();
        void SetLinearSolverBackend(LinearSolverBackend backend);
        LinearSolverBackend GetLinearSolverBackend() const;

        int GetDynamicBlockForNode(int nodeIndex) const;
        int GetDynamicBlockCount() const
        {
            return dynamicBlockCount;
        }

        const std::vector<float>& GetDeltaVelocity() const
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

    private:
        bool BuildDynamicNodeMapping(const std::vector<Node>& nodes);
        bool AssembleImplicitMatrixAndRhs(const std::vector<Node>& nodes, float dt);

        std::vector<NodeForce> nodeForces;
        std::vector<NodeMass> nodeMasses;
        std::vector<StiffnessBlock> stiffnessBlocks;

        std::vector<float> assembledMasses;
        std::vector<Vec3> assembledForces;
        std::vector<int> nodeToDynamicBlock;
        int dynamicBlockCount = 0;
        std::vector<float> rhs;
        SparseBlockMatrix matrix;
        std::vector<float> deltaVelocity;
        LinearSolverBackend linearSolverBackend = LinearSolverBackend::Current;
    };
}
