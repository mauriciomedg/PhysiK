#pragma once

#include <cstddef>
#include <map>
#include <memory>
#include <vector>

#include "PhysiK/API/Handles.h"
#include "PhysiK/Components/Component.h"
#include "PhysiK/Core/Collision/CollisionDetectionEngine.h"
#include "PhysiK/Core/Events/EventSystem.h"
#include "PhysiK/Core/PhysicsConnections/PhysicsConnection.h"
#include "PhysiK/Core/PhysicsConnections/PointConnection.h"
#include "PhysiK/Core/World/WorldState.h"
#include "PhysiK/Core/Solvers/SolverData.h"
#include "PhysiK/Core/Solvers/Linear/ConjugateGradientSolver.h"
#include "PhysiK/Math/Mat3.h"
#include "PhysiK/Math/Quaternion.h"
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
        void SetGravity(const Vec3& value);
        const Vec3& GetGravity() const;

        Node& GetNode(int index);
        const Node& GetNode(int index) const;
        const std::vector<Node>& GetNodes() const;
        Vec3& GetNodePosition(int nodeIndex);
        const Vec3& GetNodePosition(int nodeIndex) const;
        Vec3& GetNodeVelocity(int nodeIndex);
        const Vec3& GetNodeVelocity(int nodeIndex) const;
        float& GetNodeMass(int nodeIndex);
        float GetNodeMass(int nodeIndex) const;
        bool NodeHasRotation(int nodeIndex) const;
        void SetNodeHasRotation(int nodeIndex, bool hasRotation);
        Quaternion& GetNodeOrientation(int nodeIndex);
        const Quaternion& GetNodeOrientation(int nodeIndex) const;
        void SetNodeOrientation(int nodeIndex, const Quaternion& orientation);
        Vec3& GetNodeAngularVelocity(int nodeIndex);
        const Vec3& GetNodeAngularVelocity(int nodeIndex) const;
        void SetNodeAngularVelocity(int nodeIndex, const Vec3& angularVelocity);
        Vec3& GetNodeTorque(int nodeIndex);
        const Vec3& GetNodeTorque(int nodeIndex) const;
        void SetNodeTorque(int nodeIndex, const Vec3& torque);
        Mat3& GetNodeInverseInertia(int nodeIndex);
        const Mat3& GetNodeInverseInertia(int nodeIndex) const;
        void SetNodeInverseInertia(int nodeIndex, const Mat3& inverseInertia);
        void SetNodePosition(int index, const Vec3& position);
        void SetNodeFixed(int nodeIndex, bool fixed);
        bool IsNodeFixed(int nodeIndex) const;
        const std::vector<std::unique_ptr<Component>>& GetComponents() const;
        ComponentHandle GetComponentHandleByIndex(int index) const;

        int GetTransientConnectionCount() const;
        bool HasValidNodeIndices(const PointConnection& connection) const;

    private:
        bool IsComponentHandleValid(ComponentHandle handle) const;
        void RegisterComponentForExecution(Component* component);
        void UnregisterComponentFromExecution(Component* component);
        void PreUpdateComponents(float frameDt);
        void PostUpdateComponents(float frameDt);
        void BuildSolverData(SolverData& solverData, float dt);
        void ValidateNodeMasses(const SolverData& solverData) const;
        void AddGravityForces(SolverData& solverData);
        void AssembleConnectionSystems(SolverData& solverData, float dt);
        void AssembleComponentSystems(SolverData& solverData, float dt);
        void PrecomputeSolve(SolverData& solverData, float dt);
        bool SolveImplicitLinearSystem(SolverData& solverData, float dt);
        bool IntegrateImplicitEuler(const SolverData& solverData, float dt);
        void IntegrateExplicitEuler(const SolverData& solverData, float dt);
        void MarkSubstepConnectionBegin();
        void ClearSubstepConnections();
        void ClearFrameConnections();

        std::vector<Node> nodes;
        WorldState state;

        std::vector<std::unique_ptr<Component>> components;
        std::multimap<
            ComponentExecutionPriority,
            Component*,
            ComponentExecutionPriorityLess>
            orderedComponents;
        std::vector<std::uint32_t> componentGenerations;
        std::vector<std::uint32_t> freeComponentSlots;
        std::vector<std::unique_ptr<PhysicsConnection>> transientConnections;
        std::size_t firstSubstepConnectionIndex = 0u;

        CollisionDetectionEngine collisionDetectionEngine;
        EventSystem eventSystem;

        Vec3 gravity;
        int substepCount = 1;
        SolverMode solverMode = SolverMode::Explicit;
        ConjugateGradientSettings conjugateGradientSettings;
        LinearSolveResult lastConjugateGradientResult;
    };
}
