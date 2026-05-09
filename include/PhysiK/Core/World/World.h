#pragma once

#include <memory>
#include <vector>

#include "PhysiK/API/Handles.h"
#include "PhysiK/Components/Component.h"
#include "PhysiK/Core/Collision/CollisionDetectionEngine.h"
#include "PhysiK/Core/PhysicsConnections/PhysicsConnection.h"
#include "PhysiK/Core/PhysicsConnections/PointConnection.h"
#include "PhysiK/Core/Solvers/SolverData.h"
#include "PhysiK/Math/Vec3.h"
#include "PhysiK/PhysicsData/Contact.h"
#include "PhysiK/PhysicsData/Node.h"

namespace PhysiK
{
    class World
    {
    public:
        World();

        void Step(float frameDt);

        int AddNode(const Vec3& position, float inverseMass = 1.0f);
        ComponentHandle AddComponent(std::unique_ptr<Component> component);
        Component* GetComponent(ComponentHandle handle);
        const Component* GetComponent(ComponentHandle handle) const;
        void DestroyComponent(ComponentHandle handle);

        void AddPointConnection(const PointConnection& connection);
        void AddTransientConnection(std::unique_ptr<PhysicsConnection> connection);
        void SetExternalLogicCallback(ExternalLogicCallback callback, void* userData);
        void ClearExternalLogicCallback();

        void SetSubstepCount(int count);
        int GetSubstepCount() const;
        void SetGravity(const Vec3& value);
        const Vec3& GetGravity() const;

        Node& GetNode(int index);
        const Node& GetNode(int index) const;
        const std::vector<Node>& GetNodes() const;
        void SetNodePosition(int index, const Vec3& position);
        const std::vector<std::unique_ptr<Component>>& GetComponents() const;

        int GetTransientConnectionCount() const;
        bool HasValidNodeIndices(const PointConnection& connection) const;

    private:
        bool IsComponentHandleValid(ComponentHandle handle) const;
        void RunExternalLogic();
        void UpdateFrameComponents(float frameDt);
        void UpdateKinematicTargets();
        void AccumulateForces(float dt);
        void AddGravityForces(SolverData& solverData);
        void AddConnectionForces(SolverData& solverData, float dt);
        void GenerateCollisionConnections();
        void AddPhysicsModelForces(SolverData& solverData, float dt);
        void GeneratePointConnectionFromContact(const Contact& contact);
        void Solve(SolverData& solverData, float dt);
        void Integrate(float dt);
        void ClearForces();
        void ClearTransientConnections();

        std::vector<Node> nodes;

        std::vector<std::unique_ptr<Component>> components;
        std::vector<std::uint32_t> componentGenerations;
        std::vector<std::uint32_t> freeComponentSlots;
        std::vector<std::unique_ptr<PhysicsConnection>> transientConnections;

        CollisionDetectionEngine collisionDetectionEngine;

        ExternalLogicCallback externalLogicCallback = nullptr;
        void* externalLogicUserData = nullptr;

        Vec3 gravity;
        int substepCount = 1;
    };
}
