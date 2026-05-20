#pragma once

namespace PhysiK
{
    class Component;
    class World;

    enum class PhysicsEventType
    {
        TetMeshTopologyChanged
    };

    struct PhysicsEvent
    {
        PhysicsEventType type;
        World* world = nullptr;
        Component* sender = nullptr;
    };
}
