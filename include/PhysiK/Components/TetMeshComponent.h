#pragma once

#include <vector>

#include "PhysiK/Components/Component.h"
#include "PhysiK/PhysicsData/Material.h"

namespace PhysiK
{
    class TetMeshComponent : public Component
    {
    public:
        std::vector<int> nodeIndices;
        std::vector<int> tetIndices;

        Material material;
    };
}
