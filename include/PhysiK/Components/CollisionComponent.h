#pragma once

#include <vector>

#include "PhysiK/Components/Component.h"
#include "PhysiK/Math/Transform.h"

namespace PhysiK
{
    class CollisionComponent : public Component
    {
    public:
        Transform transform;

        bool generateConnections = true;
        bool generateEvents = false;
        bool isSensor = false;

        float contactStiffness = 1000.0f;
        float contactDamping = 10.0f;

        void SetKinematicTarget(const Transform& target);
        bool ConsumeKinematicTarget(Transform& outTarget);
        void PreUpdate(World& world, float dt) override;

    private:
        bool hasKinematicTarget = false;
        Transform kinematicTarget;
    };
}
