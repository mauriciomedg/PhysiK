#include "PhysiK/Core/PhysicsConnections/PointConnection.h"

#include "PhysiK/Core/Solvers/SolverData.h"
#include "PhysiK/Core/World/World.h"

#include <cmath>

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

        Mat3 ScaledIdentity(float scale)
        {
            return Mat3::FromColumns(
                Vec3{scale, 0.0f, 0.0f},
                Vec3{0.0f, scale, 0.0f},
                Vec3{0.0f, 0.0f, scale});
        }

        bool IsFinite(const Vec3& value)
        {
            return std::isfinite(value.x) &&
                std::isfinite(value.y) &&
                std::isfinite(value.z);
        }
    }

    void PointConnection::UpdateSystem(World& world, SolverData& solverData, float dt)
    {
        (void)dt;

        if (!world.HasValidNodeIndices(*this) ||
            !std::isfinite(stiffness) ||
            !std::isfinite(damping) ||
            stiffness < 0.0f ||
            damping < 0.0f)
        {
            return;
        }

        const Node& n0 = world.GetNode(node0);
        const Node& n1 = world.GetNode(node1);
        const Node& n2 = world.GetNode(node2);
        const Node& n3 = world.GetNode(node3);

        const Vec3 point = WeightedPoint(n0, n1, n2, n3, barycentric);
        const Vec3 velocity = WeightedVelocity(n0, n1, n2, n3, barycentric);
        if (!IsFinite(point) || !IsFinite(velocity) || !IsFinite(targetPosition))
        {
            return;
        }

        const Vec3 springForce = (targetPosition - point) * stiffness;
        const Vec3 dampingForce = velocity * (-damping);
        const Vec3 pointForce = springForce + dampingForce;

        const int nodeIndices[4] = {node0, node1, node2, node3};
        const float weights[4] = {
            barycentric.x,
            barycentric.y,
            barycentric.z,
            barycentric.w};

        for (int i = 0; i < 4; ++i)
        {
            if (weights[i] == 0.0f)
            {
                continue;
            }

            solverData.AddNodeForce(nodeIndices[i], pointForce * weights[i]);
        }

        if (stiffness <= 0.0f)
        {
            return;
        }

        for (int rowNode = 0; rowNode < 4; ++rowNode)
        {
            if (weights[rowNode] == 0.0f)
            {
                continue;
            }

            for (int columnNode = 0; columnNode < 4; ++columnNode)
            {
                if (weights[columnNode] == 0.0f)
                {
                    continue;
                }

                solverData.AddStiffnessBlock(
                    nodeIndices[rowNode],
                    nodeIndices[columnNode],
                    ScaledIdentity(stiffness * weights[rowNode] * weights[columnNode]));
            }
        }
    }
}
