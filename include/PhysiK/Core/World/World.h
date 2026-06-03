#pragma once

#include <memory>
#include <vector>

#include "PhysiK/API/Handles.h"
#include "PhysiK/Components/Component.h"
#include "PhysiK/Core/Collision/CollisionDetectionEngine.h"
#include "PhysiK/Core/Events/EventSystem.h"
#if defined(PHYSIK_ENABLE_PERF_LOGGING)
#include "PhysiK/Core/Performance/PerformanceLogger.h"
#endif
#include "PhysiK/Core/PhysicsConnections/PhysicsConnection.h"
#include "PhysiK/Core/PhysicsConnections/PointConnection.h"
#include "PhysiK/Core/Solvers/SolverData.h"
#include "PhysiK/Core/Solvers/Linear/ConjugateGradientSolver.h"
#include "PhysiK/Math/Vec3.h"
#include "PhysiK/PhysicsData/Node.h"

namespace PhysiK
{
    enum class SolverMode
    {
        Explicit,
        ImplicitEuler
    };

    class World
    {
    public:
        World();

        void Step(float frameDt);

        int AddNode(const Vec3& position);
        ComponentHandle AddComponent(std::unique_ptr<Component> component);
        Component* GetComponent(ComponentHandle handle);
        const Component* GetComponent(ComponentHandle handle) const;
        void DestroyComponent(ComponentHandle handle);

        void AddPointConnection(const PointConnection& connection);
        void AddTransientConnection(std::unique_ptr<PhysicsConnection> connection);
        void SubscribeToEvent(Component* listener, PhysicsEventType type);
        void UnsubscribeFromEvent(Component* listener, PhysicsEventType type);
        void EmitEvent(const PhysicsEvent& event);
        void SetExternalLogicCallback(ExternalLogicCallback callback, void* userData);
        void ClearExternalLogicCallback();

        void SetSubstepCount(int count);
        int GetSubstepCount() const;
        void SetSolverMode(SolverMode mode);
        SolverMode GetSolverMode() const;
        void SetConjugateGradientTolerance(float tolerance);
        float GetConjugateGradientTolerance() const;
        void SetConjugateGradientMaxIterations(int maxIterations);
        int GetConjugateGradientMaxIterations() const;
        int GetLastConjugateGradientIterations() const;
        float GetLastConjugateGradientResidualNorm() const;
        bool DidLastConjugateGradientSolveConverge() const;
        void EnablePerformanceLogging(bool enabled);
        void SetPerformanceLogPath(const char* path);
        void SetGravity(const Vec3& value);
        const Vec3& GetGravity() const;

        Node& GetNode(int index);
        const Node& GetNode(int index) const;
        const std::vector<Node>& GetNodes() const;
        void SetNodePosition(int index, const Vec3& position);
        void SetNodeFixed(int nodeIndex, bool fixed);
        bool IsNodeFixed(int nodeIndex) const;
        const std::vector<std::unique_ptr<Component>>& GetComponents() const;
        ComponentHandle GetComponentHandleByIndex(int index) const;

        int GetTransientConnectionCount() const;
        bool HasValidNodeIndices(const PointConnection& connection) const;

    private:
        bool IsComponentHandleValid(ComponentHandle handle) const;
        void RunExternalLogic();
        void PreUpdateComponents(float frameDt);
        void PostUpdateComponents(float frameDt);
#if defined(PHYSIK_ENABLE_PERF_LOGGING)
        void BuildSolverData(
            SolverData& solverData,
            float dt,
            PerformanceLogRecord* performanceRecord);
#else
        void BuildSolverData(SolverData& solverData, float dt);
#endif
        void ValidateNodeMasses(const SolverData& solverData) const;
        void AddGravityForces(SolverData& solverData);
        void AssembleConnectionSystems(SolverData& solverData, float dt);
        void AssembleComponentSystems(SolverData& solverData, float dt);
#if defined(PHYSIK_ENABLE_PERF_LOGGING)
        void AssembleComponentSystems(
            SolverData& solverData,
            float dt,
            PerformanceLogRecord* performanceRecord);
#endif
#if defined(PHYSIK_ENABLE_PERF_LOGGING)
        void PrecomputeSolve(
            SolverData& solverData,
            float dt,
            PerformanceLogRecord* performanceRecord);
        bool SolveImplicitLinearSystem(
            SolverData& solverData,
            float dt,
            PerformanceLogRecord* performanceRecord);
#else
        void PrecomputeSolve(SolverData& solverData, float dt);
        bool SolveImplicitLinearSystem(SolverData& solverData, float dt);
#endif
        bool IntegrateImplicitEuler(const SolverData& solverData, float dt);
        void IntegrateExplicitEuler(const SolverData& solverData, float dt);
        void ClearTransientConnections();

        std::vector<Node> nodes;

        std::vector<std::unique_ptr<Component>> components;
        std::vector<std::uint32_t> componentGenerations;
        std::vector<std::uint32_t> freeComponentSlots;
        std::vector<std::unique_ptr<PhysicsConnection>> transientConnections;

        CollisionDetectionEngine collisionDetectionEngine;
        EventSystem eventSystem;

        ExternalLogicCallback externalLogicCallback = nullptr;
        void* externalLogicUserData = nullptr;

        Vec3 gravity;
        int substepCount = 1;
        SolverMode solverMode = SolverMode::Explicit;
        ConjugateGradientSettings conjugateGradientSettings;
        ConjugateGradientResult lastConjugateGradientResult;
#if defined(PHYSIK_ENABLE_PERF_LOGGING)
        std::uint64_t frameIndex = 0;
        PerformanceLogger performanceLogger;
#endif
    };
}
