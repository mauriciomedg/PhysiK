#include "PhysiK/Core/Physics/FEM/FEMModel.h"

#include <algorithm>
#include <cmath>

namespace PhysiK
{
    namespace
    {
        Mat3 BuildDm(const Tet& tet, const std::vector<Node>& nodes)
        {
            const Vec3& x0 = nodes[static_cast<std::size_t>(tet.node0)].position;
            const Vec3& x1 = nodes[static_cast<std::size_t>(tet.node1)].position;
            const Vec3& x2 = nodes[static_cast<std::size_t>(tet.node2)].position;
            const Vec3& x3 = nodes[static_cast<std::size_t>(tet.node3)].position;
            return Mat3::FromColumns(x1 - x0, x2 - x0, x3 - x0);
        }

        bool HasValidNodes(const Tet& tet, const std::vector<Node>& nodes)
        {
            const int nodeCount = static_cast<int>(nodes.size());
            return tet.node0 >= 0 && tet.node0 < nodeCount &&
                tet.node1 >= 0 && tet.node1 < nodeCount &&
                tet.node2 >= 0 && tet.node2 < nodeCount &&
                tet.node3 >= 0 && tet.node3 < nodeCount;
        }
    }

    void FEMModel::InitializeTetRestData(Tet& tet, const std::vector<Node>& nodes)
    {
        if (!HasValidNodes(tet, nodes))
        {
            return;
        }

        const Mat3 restDm = BuildDm(tet, nodes);
        const float determinant = Determinant(restDm);
        tet.restVolume = std::abs(determinant) / 6.0f;
        tet.restDmInverse = Inverse(restDm);
    }

    void FEMModel::AccumulateElasticForces(const std::vector<Tet>& tets, std::vector<Node>& nodes)
    {
        for (const Tet& tet : tets)
        {
            if (!HasValidNodes(tet, nodes) || tet.restVolume <= 0.0f)
            {
                continue;
            }

            const Mat3 deformationGradient = BuildDm(tet, nodes) * tet.restDmInverse;
            const Mat3 strain = deformationGradient - Mat3::Identity();
            const Mat3 forceMatrix =
                strain * Transpose(tet.restDmInverse);

            const Vec3 force1 = forceMatrix.columns[0] * (-tet.stiffness * tet.restVolume);
            const Vec3 force2 = forceMatrix.columns[1] * (-tet.stiffness * tet.restVolume);
            const Vec3 force3 = forceMatrix.columns[2] * (-tet.stiffness * tet.restVolume);
            const Vec3 force0 = -(force1 + force2 + force3);

            Node& node0 = nodes[static_cast<std::size_t>(tet.node0)];
            Node& node1 = nodes[static_cast<std::size_t>(tet.node1)];
            Node& node2 = nodes[static_cast<std::size_t>(tet.node2)];
            Node& node3 = nodes[static_cast<std::size_t>(tet.node3)];

            node0.force += force0 - node0.velocity * tet.damping;
            node1.force += force1 - node1.velocity * tet.damping;
            node2.force += force2 - node2.velocity * tet.damping;
            node3.force += force3 - node3.velocity * tet.damping;
        }
    }
}
