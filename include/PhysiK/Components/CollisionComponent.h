#pragma once

#include <vector>

#include "PhysiK/Components/Component.h"
#include "PhysiK/Math/Transform.h"
#include "PhysiK/PhysicsData/Contact.h"

namespace PhysiK
{
    class CollisionDetectionEngine;

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

        virtual void QueryContacts(
            World& world,
            CollisionDetectionEngine& collisionDetectionEngine,
            std::vector<Contact>& outContacts) = 0;

    private:
        bool hasKinematicTarget = false;
        Transform kinematicTarget;
    };
}
