#pragma once

namespace PhysiK
{
    struct Material
    {
        float density = 1.0f;
        float youngModulus = 25.0f;
        float poissonRatio = 0.3f;
        float damping = 0.25f;
    };
}
