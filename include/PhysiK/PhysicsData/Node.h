#pragma once

namespace PhysiK
{
    struct Node
    {
        int stateIndex = -1;
        bool active = true;
        bool fixed = false;
        bool hasRotation = false;
    };
}
