#include "PhysiK/Components/Component.h"
#include "PhysiK/Components/ComponentExecutionPriority.h"

#include <cstdio>
#include <cstdlib>
#include <map>
#include <vector>

namespace
{
    class DefaultPriorityTestComponent :
        public PhysiK::Component
    {
    };

    void Require(bool condition, const char* message)
    {
        if (!condition)
        {
            std::fprintf(stderr, "ComponentExecutionPriorityTests failed: %s\n", message);
            std::exit(1);
        }
    }
}

void DefaultComponentPriorityIsDefault()
{
    DefaultPriorityTestComponent component;

    Require(
        component.GetExecutionPriority() ==
            PhysiK::ComponentExecutionPriority::
                Default,
        "default component priority should be Default");
}

void ComponentExecutionPriorityLessOrdersValues()
{
    const PhysiK::ComponentExecutionPriorityLess less;

    Require(
        less(
            PhysiK::ComponentExecutionPriority::
                ScriptComponent,
            PhysiK::ComponentExecutionPriority::
                CollisionSphereComponent),
        "ScriptComponent should sort before CollisionSphereComponent");
    Require(
        !less(
            PhysiK::ComponentExecutionPriority::
                VisualMeshComponent,
            PhysiK::ComponentExecutionPriority::
                ScriptComponent),
        "VisualMeshComponent should not sort before ScriptComponent");
}

void ComponentExecutionPriorityMultimapPreservesPriorityOrderAndDuplicates()
{
    using PhysiK::ComponentExecutionPriority;

    std::multimap<
        ComponentExecutionPriority,
        int,
        PhysiK::ComponentExecutionPriorityLess>
        orderedValues;

    orderedValues.insert({ComponentExecutionPriority::Default, 10});
    orderedValues.insert({ComponentExecutionPriority::TetMeshPhysicsComponent, 4});
    orderedValues.insert({ComponentExecutionPriority::CollisionSphereComponent, 2});
    orderedValues.insert({ComponentExecutionPriority::VisualMeshComponent, 9});
    orderedValues.insert({ComponentExecutionPriority::ScriptComponent, 1});
    orderedValues.insert({ComponentExecutionPriority::SurfaceVisualComponent, 8});
    orderedValues.insert({ComponentExecutionPriority::TetMeshComponent, 3});
    orderedValues.insert({ComponentExecutionPriority::TetMeshMapperComponent, 6});
    orderedValues.insert({ComponentExecutionPriority::TopologyMeshComponent, 5});
    orderedValues.insert({ComponentExecutionPriority::SurfaceExtractionComponent, 7});
    orderedValues.insert({ComponentExecutionPriority::CollisionSphereComponent, 20});

    std::vector<ComponentExecutionPriority> priorities;
    std::vector<int> collisionSphereValues;
    for (const auto& entry : orderedValues)
    {
        priorities.push_back(entry.first);
        if (entry.first == ComponentExecutionPriority::CollisionSphereComponent)
        {
            collisionSphereValues.push_back(entry.second);
        }
    }

    const std::vector<ComponentExecutionPriority> expectedPriorities{
        ComponentExecutionPriority::ScriptComponent,
        ComponentExecutionPriority::CollisionSphereComponent,
        ComponentExecutionPriority::CollisionSphereComponent,
        ComponentExecutionPriority::TetMeshComponent,
        ComponentExecutionPriority::TetMeshPhysicsComponent,
        ComponentExecutionPriority::TopologyMeshComponent,
        ComponentExecutionPriority::TetMeshMapperComponent,
        ComponentExecutionPriority::SurfaceExtractionComponent,
        ComponentExecutionPriority::SurfaceVisualComponent,
        ComponentExecutionPriority::VisualMeshComponent,
        ComponentExecutionPriority::Default};

    Require(
        priorities == expectedPriorities,
        "multimap iteration should follow component execution priority order");
    Require(
        collisionSphereValues.size() == 2u,
        "multimap should preserve duplicate CollisionSphereComponent entries");
    Require(
        collisionSphereValues[0] == 2 &&
            collisionSphereValues[1] == 20,
        "duplicate CollisionSphereComponent entries should preserve inserted values");
}

int main()
{
    DefaultComponentPriorityIsDefault();
    ComponentExecutionPriorityLessOrdersValues();
    ComponentExecutionPriorityMultimapPreservesPriorityOrderAndDuplicates();

    return 0;
}
