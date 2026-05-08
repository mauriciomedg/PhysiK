#include "PhysiK/Core/PhysicsConnections/PointConnection.h"

#include "PhysiK/Core/Solvers/SolverData.h"
#include "PhysiK/Core/World/World.h"

namespace PhysiK
{
    namespace
    {
        Vec3 WeightedPoint(
            const Node& node0,
            const Node& node1,
            const Node& node2,
            const Node& node3,
            const Vec4& weights)
        {
            return node0.position * weights.x +
                node1.position * weights.y +
                node2.position * weights.z +
                node3.position * weights.w;
        }

        Vec3 WeightedVelocity(
            const Node& node0,
            const Node& node1,
            const Node& node2,
            const Node& node3,
            const Vec4& weights)
        {
            return node0.velocity * weights.x +
                node1.velocity * weights.y +
                node2.velocity * weights.z +
                node3.velocity * weights.w;
        }
    }

    void PointConnection::UpdateSystem(World& world, SolverData& solverData, float dt)
    {
        (void)dt;

        if (!world.HasValidNodeIndices(*this))
        {
            return;
        }

        const Node& n0 = world.GetNode(node0);
        const Node& n1 = world.GetNode(node1);
        const Node& n2 = world.GetNode(node2);
        const Node& n3 = world.GetNode(node3);

        const Vec3 point = WeightedPoint(n0, n1, n2, n3, barycentric);
        const Vec3 velocity = WeightedVelocity(n0, n1, n2, n3, barycentric);
        const Vec3 springForce = (targetPosition - point) * stiffness;
        const Vec3 dampingForce = velocity * (-damping);
        const Vec3 pointForce = springForce + dampingForce;

        solverData.AddNodeForce(node0, pointForce * barycentric.x);
        solverData.AddNodeForce(node1, pointForce * barycentric.y);
        solverData.AddNodeForce(node2, pointForce * barycentric.z);
        solverData.AddNodeForce(node3, pointForce * barycentric.w);
    }
}
