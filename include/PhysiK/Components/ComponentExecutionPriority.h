#pragma once

namespace PhysiK
{
    enum class ComponentExecutionPriority : int
    {
        ScriptComponent =
            100,

        CollisionSphereComponent =
            200,

        TetMeshComponent =
            300,

        TetMeshPhysicsComponent =
            400,

        TopologyMeshComponent =
            500,

        TetMeshMapperComponent =
            600,

        SurfaceExtractionComponent =
            700,

        SurfaceVisualComponent =
            800,

        VisualMeshComponent =
            900,

        Default =
            1000
    };

    struct ComponentExecutionPriorityLess
    {
        bool operator()(
            ComponentExecutionPriority left,
            ComponentExecutionPriority right) const
        {
            return static_cast<int>(
                       left) <
                static_cast<int>(
                       right);
        }
    };
}
