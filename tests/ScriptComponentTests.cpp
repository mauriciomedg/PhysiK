#include "PhysiK/API/PhysiKAPI.h"
#include "PhysiK/Components/ScriptComponent.h"
#include "PhysiK/Core/World/World.h"

#include <cstdio>
#include <cstdlib>

namespace
{
    struct ScriptCallbackTestContext
    {
        int callCount = 0;
        PhysiK::WorldHandle receivedWorld = nullptr;
        void* receivedUserData = nullptr;
    };

    void Require(bool condition, const char* message)
    {
        if (!condition)
        {
            std::fprintf(stderr, "ScriptComponentTests failed: %s\n", message);
            std::exit(1);
        }
    }

    void RecordScriptCallback(
        PhysiK::WorldHandle world,
        void* userData)
    {
        auto* context =
            static_cast<ScriptCallbackTestContext*>(userData);

        ++context->callCount;
        context->receivedWorld =
            world;
        context->receivedUserData =
            userData;
    }

    PhysiK::ScriptComponent* GetScriptComponent(
        PhysiK::WorldHandle worldHandle,
        PhysiK::ComponentHandle componentHandle)
    {
        auto* world =
            static_cast<PhysiK::World*>(worldHandle);
        if (world == nullptr)
        {
            return nullptr;
        }

        return dynamic_cast<PhysiK::ScriptComponent*>(
            world->GetComponent(componentHandle));
    }
}

void ScriptComponentDefaultStateHasNoCallback()
{
    PhysiK::ScriptComponent component;

    Require(
        component.GetExternalLogicCallback() == nullptr,
        "default callback should be null");
    Require(
        component.GetExternalLogicUserData() == nullptr,
        "default user data should be null");
}

void ScriptComponentPreUpdateWithoutCallbackDoesNothing()
{
    PhysiK::WorldHandle worldHandle =
        PHYSIK_CreateWorld();
    Require(
        worldHandle != nullptr,
        "world creation failed for no-callback test");

    auto& world =
        *static_cast<PhysiK::World*>(worldHandle);
    PhysiK::ScriptComponent component;

    component.PreUpdate(
        world,
        1.0f / 60.0f);

    PHYSIK_DestroyWorld(
        worldHandle);
}

void ScriptComponentPreUpdateInvokesCallback()
{
    PhysiK::WorldHandle worldHandle =
        PHYSIK_CreateWorld();
    Require(
        worldHandle != nullptr,
        "world creation failed for callback test");

    auto& world =
        *static_cast<PhysiK::World*>(worldHandle);
    PhysiK::ScriptComponent component;
    ScriptCallbackTestContext context;

    component.SetExternalLogicCallback(
        RecordScriptCallback,
        &context);

    component.PreUpdate(
        world,
        1.0f / 60.0f);

    Require(
        context.callCount == 1,
        "callback should execute exactly once");
    Require(
        context.receivedWorld == static_cast<PhysiK::WorldHandle>(&world),
        "callback should receive the same world address");
    Require(
        context.receivedUserData == &context,
        "callback should receive the configured user data");

    PHYSIK_DestroyWorld(
        worldHandle);
}

void ScriptComponentCallbackRunsOncePerPreUpdate()
{
    PhysiK::WorldHandle worldHandle =
        PHYSIK_CreateWorld();
    Require(
        worldHandle != nullptr,
        "world creation failed for repeated-callback test");

    auto& world =
        *static_cast<PhysiK::World*>(worldHandle);
    PhysiK::ScriptComponent component;
    ScriptCallbackTestContext context;

    component.SetExternalLogicCallback(
        RecordScriptCallback,
        &context);

    component.PreUpdate(
        world,
        1.0f / 60.0f);
    component.PreUpdate(
        world,
        1.0f / 60.0f);
    component.PreUpdate(
        world,
        1.0f / 60.0f);

    Require(
        context.callCount == 3,
        "callback should execute once per PreUpdate call");

    PHYSIK_DestroyWorld(
        worldHandle);
}

void ScriptComponentClearCallbackPreventsInvocation()
{
    PhysiK::WorldHandle worldHandle =
        PHYSIK_CreateWorld();
    Require(
        worldHandle != nullptr,
        "world creation failed for clear-callback test");

    auto& world =
        *static_cast<PhysiK::World*>(worldHandle);
    PhysiK::ScriptComponent component;
    ScriptCallbackTestContext context;

    component.SetExternalLogicCallback(
        RecordScriptCallback,
        &context);
    component.ClearExternalLogicCallback();

    component.PreUpdate(
        world,
        1.0f / 60.0f);

    Require(
        context.callCount == 0,
        "cleared callback should not execute");
    Require(
        component.GetExternalLogicCallback() == nullptr,
        "cleared callback should be null");
    Require(
        component.GetExternalLogicUserData() == nullptr,
        "cleared user data should be null");

    PHYSIK_DestroyWorld(
        worldHandle);
}

void ScriptComponentReplacingCallbackUserDataUsesNewContext()
{
    PhysiK::WorldHandle worldHandle =
        PHYSIK_CreateWorld();
    Require(
        worldHandle != nullptr,
        "world creation failed for replace-callback test");

    auto& world =
        *static_cast<PhysiK::World*>(worldHandle);
    PhysiK::ScriptComponent component;
    ScriptCallbackTestContext firstContext;
    ScriptCallbackTestContext secondContext;

    component.SetExternalLogicCallback(
        RecordScriptCallback,
        &firstContext);
    component.PreUpdate(
        world,
        1.0f / 60.0f);

    component.SetExternalLogicCallback(
        RecordScriptCallback,
        &secondContext);
    component.PreUpdate(
        world,
        1.0f / 60.0f);

    Require(
        firstContext.callCount == 1,
        "first callback context should execute once");
    Require(
        secondContext.callCount == 1,
        "second callback context should execute once");
    Require(
        secondContext.receivedUserData == &secondContext,
        "replacement callback should receive the second user data");

    PHYSIK_DestroyWorld(
        worldHandle);
}

void ScriptComponentApiCreatesValidComponent()
{
    PhysiK::WorldHandle world =
        PHYSIK_CreateWorld();
    Require(
        world != nullptr,
        "world creation failed for script API create test");

    const PhysiK::ComponentHandle script =
        PHYSIK_CreateScriptComponent(world);

    Require(
        PHYSIK_IsComponentHandleValid(world, script) == 1,
        "created ScriptComponent handle should be valid");

    PHYSIK_DestroyWorld(world);
}

void ScriptComponentApiConfiguresCallback()
{
    PhysiK::WorldHandle world =
        PHYSIK_CreateWorld();
    Require(
        world != nullptr,
        "world creation failed for script API callback test");

    const PhysiK::ComponentHandle script =
        PHYSIK_CreateScriptComponent(world);
    ScriptCallbackTestContext context;

    PHYSIK_SetScriptComponentCallback(
        world,
        script,
        RecordScriptCallback,
        &context);

    PhysiK::ScriptComponent* component =
        GetScriptComponent(world, script);
    Require(
        component != nullptr,
        "created script component should be retrievable");

    component->PreUpdate(
        *static_cast<PhysiK::World*>(world),
        1.0f / 60.0f);

    Require(
        context.callCount == 1,
        "script API callback should execute exactly once");
    Require(
        context.receivedWorld == world,
        "script API callback should receive the world handle");
    Require(
        context.receivedUserData == &context,
        "script API callback should receive user data");

    PHYSIK_DestroyWorld(world);
}

void ScriptComponentApiClearsCallback()
{
    PhysiK::WorldHandle world =
        PHYSIK_CreateWorld();
    Require(
        world != nullptr,
        "world creation failed for script API clear test");

    const PhysiK::ComponentHandle script =
        PHYSIK_CreateScriptComponent(world);
    ScriptCallbackTestContext context;

    PHYSIK_SetScriptComponentCallback(
        world,
        script,
        RecordScriptCallback,
        &context);
    PHYSIK_ClearScriptComponentCallback(
        world,
        script);

    PhysiK::ScriptComponent* component =
        GetScriptComponent(world, script);
    Require(
        component != nullptr,
        "script component should be retrievable after clearing callback");

    component->PreUpdate(
        *static_cast<PhysiK::World*>(world),
        1.0f / 60.0f);

    Require(
        context.callCount == 0,
        "cleared script API callback should not execute");

    PHYSIK_DestroyWorld(world);
}

void ScriptComponentApiSupportsMultipleIndependentComponents()
{
    PhysiK::WorldHandle world =
        PHYSIK_CreateWorld();
    Require(
        world != nullptr,
        "world creation failed for multiple script API test");

    const PhysiK::ComponentHandle scriptA =
        PHYSIK_CreateScriptComponent(world);
    const PhysiK::ComponentHandle scriptB =
        PHYSIK_CreateScriptComponent(world);
    const PhysiK::ComponentHandle scriptC =
        PHYSIK_CreateScriptComponent(world);

    ScriptCallbackTestContext contextA;
    ScriptCallbackTestContext contextB;
    ScriptCallbackTestContext contextC;

    PHYSIK_SetScriptComponentCallback(
        world,
        scriptA,
        RecordScriptCallback,
        &contextA);
    PHYSIK_SetScriptComponentCallback(
        world,
        scriptB,
        RecordScriptCallback,
        &contextB);
    PHYSIK_SetScriptComponentCallback(
        world,
        scriptC,
        RecordScriptCallback,
        &contextC);

    PhysiK::ScriptComponent* componentA =
        GetScriptComponent(world, scriptA);
    PhysiK::ScriptComponent* componentB =
        GetScriptComponent(world, scriptB);
    PhysiK::ScriptComponent* componentC =
        GetScriptComponent(world, scriptC);
    Require(
        componentA != nullptr &&
            componentB != nullptr &&
            componentC != nullptr,
        "all script components should be retrievable");

    auto& resolvedWorld =
        *static_cast<PhysiK::World*>(world);
    componentA->PreUpdate(
        resolvedWorld,
        1.0f / 60.0f);
    componentB->PreUpdate(
        resolvedWorld,
        1.0f / 60.0f);
    componentC->PreUpdate(
        resolvedWorld,
        1.0f / 60.0f);

    Require(
        contextA.callCount == 1,
        "script component A callback should execute once");
    Require(
        contextB.callCount == 1,
        "script component B callback should execute once");
    Require(
        contextC.callCount == 1,
        "script component C callback should execute once");

    PHYSIK_DestroyWorld(world);
}

void ScriptComponentApiUsesGenericDestroyComponent()
{
    PhysiK::WorldHandle world =
        PHYSIK_CreateWorld();
    Require(
        world != nullptr,
        "world creation failed for script API destroy test");

    const PhysiK::ComponentHandle scriptA =
        PHYSIK_CreateScriptComponent(world);
    const PhysiK::ComponentHandle scriptB =
        PHYSIK_CreateScriptComponent(world);

    PHYSIK_DestroyComponent(
        world,
        scriptA);

    Require(
        PHYSIK_IsComponentHandleValid(world, scriptA) == 0,
        "destroyed script component handle should be invalid");
    Require(
        PHYSIK_IsComponentHandleValid(world, scriptB) == 1,
        "other script component handle should remain valid");

    PHYSIK_DestroyWorld(world);
}

void ScriptComponentApiHandlesInvalidScriptHandle()
{
    PhysiK::WorldHandle world =
        PHYSIK_CreateWorld();
    Require(
        world != nullptr,
        "world creation failed for invalid script API test");

    const PhysiK::ComponentHandle script =
        PHYSIK_CreateScriptComponent(world);
    ScriptCallbackTestContext validContext;
    ScriptCallbackTestContext invalidContext;

    PHYSIK_SetScriptComponentCallback(
        world,
        script,
        RecordScriptCallback,
        &validContext);

    PHYSIK_SetScriptComponentCallback(
        world,
        PhysiK::ComponentHandle{},
        RecordScriptCallback,
        &invalidContext);
    PHYSIK_ClearScriptComponentCallback(
        world,
        PhysiK::ComponentHandle{});

    PhysiK::ScriptComponent* component =
        GetScriptComponent(world, script);
    Require(
        component != nullptr,
        "valid script component should remain retrievable after invalid-handle calls");

    component->PreUpdate(
        *static_cast<PhysiK::World*>(world),
        1.0f / 60.0f);

    Require(
        validContext.callCount == 1,
        "invalid handle calls should not modify unrelated script component");
    Require(
        invalidContext.callCount == 0,
        "invalid handle callback should not execute");

    PHYSIK_DestroyWorld(world);
}

void ScriptComponentApiIgnoresWrongComponentType()
{
    PhysiK::WorldHandle world =
        PHYSIK_CreateWorld();
    Require(
        world != nullptr,
        "world creation failed for wrong-type script API test");

    const PhysiK::ComponentHandle script =
        PHYSIK_CreateScriptComponent(world);
    const PhysiK::ComponentHandle sphere =
        PHYSIK_CreateCollisionSphereComponent(
            world,
            0.0f,
            0.0f,
            0.0f,
            1.0f);

    ScriptCallbackTestContext scriptContext;
    ScriptCallbackTestContext wrongTypeContext;

    PHYSIK_SetScriptComponentCallback(
        world,
        script,
        RecordScriptCallback,
        &scriptContext);
    PHYSIK_SetScriptComponentCallback(
        world,
        sphere,
        RecordScriptCallback,
        &wrongTypeContext);
    PHYSIK_ClearScriptComponentCallback(
        world,
        sphere);

    PhysiK::ScriptComponent* component =
        GetScriptComponent(world, script);
    Require(
        component != nullptr,
        "script component should remain retrievable after wrong-type calls");

    component->PreUpdate(
        *static_cast<PhysiK::World*>(world),
        1.0f / 60.0f);

    Require(
        scriptContext.callCount == 1,
        "wrong-type calls should not modify unrelated script component");
    Require(
        wrongTypeContext.callCount == 0,
        "wrong-type callback should not execute");

    PHYSIK_DestroyWorld(world);
}

int main()
{
    ScriptComponentDefaultStateHasNoCallback();
    ScriptComponentPreUpdateWithoutCallbackDoesNothing();
    ScriptComponentPreUpdateInvokesCallback();
    ScriptComponentCallbackRunsOncePerPreUpdate();
    ScriptComponentClearCallbackPreventsInvocation();
    ScriptComponentReplacingCallbackUserDataUsesNewContext();
    ScriptComponentApiCreatesValidComponent();
    ScriptComponentApiConfiguresCallback();
    ScriptComponentApiClearsCallback();
    ScriptComponentApiSupportsMultipleIndependentComponents();
    ScriptComponentApiUsesGenericDestroyComponent();
    ScriptComponentApiHandlesInvalidScriptHandle();
    ScriptComponentApiIgnoresWrongComponentType();

    return 0;
}
