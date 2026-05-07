#pragma once

#include "PhysiK/Components/CollisionComponent.h"

namespace PhysiK
{
    class CollisionSphereComponent : public CollisionComponent
    {
    public:
        float radius = 0.5f;

        void QueryContacts(World& world, std::vector<Contact>& outContacts) override;
    };
}
