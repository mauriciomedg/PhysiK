#pragma once

#include <cstdint>

namespace PhysiK
{
    using WorldHandle = void*;
    using GeneratedTetMeshHandle = void*;

    struct ComponentHandle
    {
        std::uint32_t index = 0xFFFFFFFFu;
        std::uint32_t generation = 0u;

        bool IsValid() const
        {
            return index != 0xFFFFFFFFu && generation != 0u;
        }
    };

    using RigidBodyHandle = ComponentHandle;

    using ExternalLogicCallback = void (*)(WorldHandle world, void* userData);
}
