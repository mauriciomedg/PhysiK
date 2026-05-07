#include "PhysiK/Core/Collision/CollisionDetectionEngine.h"

#include "PhysiK/Components/CollisionComponent.h"

namespace PhysiK
{
    void CollisionDetectionEngine::QueryContacts(
        World& world,
        CollisionComponent& component,
        std::vector<Contact>& outContacts)
    {
        component.QueryContacts(world, outContacts);
    }
}
