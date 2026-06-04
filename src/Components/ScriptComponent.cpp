#include "PhysiK/Components/ScriptComponent.h"

#include "PhysiK/Core/World/World.h"

namespace PhysiK
{
    void ScriptComponent::SetExternalLogicCallback(
        ExternalLogicCallback callback,
        void* userData)
    {
        externalLogicCallback =
            callback;

        externalLogicUserData =
            userData;
    }

    void ScriptComponent::ClearExternalLogicCallback()
    {
        externalLogicCallback =
            nullptr;

        externalLogicUserData =
            nullptr;
    }

    ExternalLogicCallback
    ScriptComponent::GetExternalLogicCallback() const
    {
        return externalLogicCallback;
    }

    void*
    ScriptComponent::GetExternalLogicUserData() const
    {
        return externalLogicUserData;
    }

    void ScriptComponent::PreUpdate(
        World& world,
        float frameDt)
    {
        (void)frameDt;

        if (externalLogicCallback ==
            nullptr)
        {
            return;
        }

        externalLogicCallback(
            static_cast<WorldHandle>(
                &world),
            externalLogicUserData);
    }
}
