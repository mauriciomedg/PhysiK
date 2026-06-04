#include "PhysiK/API/PhysiKAPI.h"
#include "PhysiK/Components/ScriptComponent.h"
#include "PhysiK/Core/World/World.h"

#include <cstdio>
#include <cstdlib>

namespace
{
    struct CallbackState
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

    void CountCallback(
        PhysiK::WorldHandle world,
        void* userData)
    {
        auto* state =
            static_cast<CallbackState*>(userData);

        ++state->callCount;
        state->receivedWorld =
            world;
        state->receivedUserData =
            userData;
    }
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
        0.016f);

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
    CallbackState state;

    component.SetExternalLogicCallback(
        CountCallback,
        &state);

    component.PreUpdate(
        world,
        0.016f);

    Require(
        state.callCount == 1,
        "callback should execute exactly once");
    Require(
        state.receivedWorld == static_cast<PhysiK::WorldHandle>(&world),
        "callback should receive the same world address");
    Require(
        state.receivedUserData == &state,
        "callback should receive the configured user data");

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
    CallbackState state;

    component.SetExternalLogicCallback(
        CountCallback,
        &state);
    component.ClearExternalLogicCallback();

    component.PreUpdate(
        world,
        0.016f);

    Require(
        state.callCount == 0,
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

int main()
{
    ScriptComponentPreUpdateWithoutCallbackDoesNothing();
    ScriptComponentPreUpdateInvokesCallback();
    ScriptComponentClearCallbackPreventsInvocation();

    return 0;
}
