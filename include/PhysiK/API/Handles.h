#pragma once

namespace PhysiK
{
    using WorldHandle = void*;
    using ComponentHandle = void*;
    using ExternalLogicCallback = void (*)(WorldHandle world, void* userData);
}
