#pragma once

#include <vector>

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

        void Clear()
        {
            nodeForces.clear();
        }

        void AddNodeForce(int node, const Vec3& force)
        {
            nodeForces.push_back(NodeForce{node, force});
        }

        const std::vector<NodeForce>& GetNodeForces() const
        {
            return nodeForces;
        }

    private:
        std::vector<NodeForce> nodeForces;
    };
}
