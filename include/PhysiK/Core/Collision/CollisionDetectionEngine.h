#pragma once

#include <vector>

#include "PhysiK/PhysicsData/Contact.h"

namespace PhysiK
{
    class CollisionComponent;
    class World;

    class CollisionDetectionEngine
    {
    public:
        void QueryContacts(
            World& world,
            CollisionComponent& component,
            std::vector<Contact>& outContacts);
    };
}
