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

int main()
{
    ScriptComponentDefaultStateHasNoCallback();
    ScriptComponentPreUpdateWithoutCallbackDoesNothing();
    ScriptComponentPreUpdateInvokesCallback();
    ScriptComponentCallbackRunsOncePerPreUpdate();
    ScriptComponentClearCallbackPreventsInvocation();
    ScriptComponentReplacingCallbackUserDataUsesNewContext();

    return 0;
}
