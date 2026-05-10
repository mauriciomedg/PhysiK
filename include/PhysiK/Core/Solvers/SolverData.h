#pragma once

#include <vector>

#include "PhysiK/Math/Mat3.h"
#include "PhysiK/Math/Vec3.h"

namespace PhysiK
{
    class SolverData
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

        void Clear()
        {
            nodeForces.clear();
            nodeMasses.clear();
            stiffnessBlocks.clear();
        }

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

    private:
        std::vector<NodeForce> nodeForces;
        std::vector<NodeMass> nodeMasses;
        std::vector<StiffnessBlock> stiffnessBlocks;
    };
}
