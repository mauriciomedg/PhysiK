#pragma once

#include "PhysiK/API/Handles.h"
#include "PhysiK/Components/Component.h"

namespace PhysiK
{
    class ScriptComponent final :
        public Component
    {
    public:
        void SetExternalLogicCallback(
            ExternalLogicCallback callback,
            void* userData);

        void ClearExternalLogicCallback();

        ExternalLogicCallback
        GetExternalLogicCallback() const;

        void*
        GetExternalLogicUserData() const;

        void PreUpdate(
            World& world,
            float frameDt) override;

    private:
        ExternalLogicCallback externalLogicCallback =
            nullptr;

        void* externalLogicUserData =
            nullptr;
    };
}
