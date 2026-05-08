#pragma once

#include <vector>

#include "PhysiK/PhysicsData/Node.h"
#include "PhysiK/PhysicsData/Tet.h"

namespace PhysiK
{
    class FEMModel
    {
    public:
        static void InitializeTetRestData(Tet& tet, const std::vector<Node>& nodes);
        static void AccumulateElasticForces(const std::vector<Tet>& tets, std::vector<Node>& nodes);
    };
}
