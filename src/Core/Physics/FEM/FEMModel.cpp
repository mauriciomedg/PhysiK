#include "PhysiK/Core/Physics/FEM/FEMModel.h"

#include <algorithm>
#include <array>
#include <cmath>

namespace PhysiK
{
    namespace
    {
        using Vector6 = std::array<float, 6>;
        using Vector12 = std::array<float, 12>;

        constexpr float MinTetVolume = 1.0e-8f;
        constexpr float MinPolarDeterminant = 1.0e-8f;

        Mat3 BuildDm(const Tet& tet, const std::vector<Vec3>& positions)
        {
            const Vec3& x0 = positions[static_cast<std::size_t>(tet.node0)];
            const Vec3& x1 = positions[static_cast<std::size_t>(tet.node1)];
            const Vec3& x2 = positions[static_cast<std::size_t>(tet.node2)];
            const Vec3& x3 = positions[static_cast<std::size_t>(tet.node3)];
            return Mat3::FromColumns(x1 - x0, x2 - x0, x3 - x0);
        }

        Vec3 GetRow(const Mat3& matrix, int row)
        {
            if (row == 0)
            {
                return Vec3{matrix.columns[0].x, matrix.columns[1].x, matrix.columns[2].x};
            }

            if (row == 1)
            {
                return Vec3{matrix.columns[0].y, matrix.columns[1].y, matrix.columns[2].y};
            }

            return Vec3{matrix.columns[0].z, matrix.columns[1].z, matrix.columns[2].z};
        }

        float GetBlockValue(const Mat3& matrix, int row, int column)
        {
            const Vec3& sourceColumn = matrix.columns[column];
            if (row == 0)
            {
                return sourceColumn.x;
            }

            if (row == 1)
            {
                return sourceColumn.y;
            }

            return sourceColumn.z;
        }

        void SetBlockValue(Mat3& matrix, int row, int column, float value)
        {
            Vec3& targetColumn = matrix.columns[column];
            if (row == 0)
            {
                targetColumn.x = value;
            }
            else if (row == 1)
            {
                targetColumn.y = value;
            }
            else
            {
                targetColumn.z = value;
            }
        }

        Mat3 Add(const Mat3& a, const Mat3& b)
        {
            return Mat3::FromColumns(
                a.columns[0] + b.columns[0],
                a.columns[1] + b.columns[1],
                a.columns[2] + b.columns[2]);
        }

        Mat3 Scale(const Mat3& matrix, float scale)
        {
            return Mat3::FromColumns(
                matrix.columns[0] * scale,
                matrix.columns[1] * scale,
                matrix.columns[2] * scale);
        }

        bool IsFinite(const Vec3& value)
        {
            return std::isfinite(value.x) &&
                std::isfinite(value.y) &&
                std::isfinite(value.z);
        }

        bool IsFinite(const Mat3& matrix)
        {
            return IsFinite(matrix.columns[0]) &&
                IsFinite(matrix.columns[1]) &&
                IsFinite(matrix.columns[2]);
        }

        float FrobeniusDifferenceSquared(const Mat3& a, const Mat3& b)
        {
            float total = 0.0f;
            for (int column = 0; column < 3; ++column)
            {
                const Vec3 difference = a.columns[column] - b.columns[column];
                total += difference.LengthSquared();
            }
            return total;
        }

        Mat3 ExtractRotationPolar(const Mat3& deformationGradient)
        {
            if (!IsFinite(deformationGradient) ||
                std::abs(Determinant(deformationGradient)) <= MinPolarDeterminant)
            {
                return Mat3::Identity();
            }

            Mat3 rotation = deformationGradient;
            int iterationCount = 0;
            bool earlyExit = false;
            for (int iteration = 0; iteration < 12; ++iteration)
            {
                iterationCount = iteration + 1;
                const Mat3 transposeInverse = Inverse(Transpose(rotation));
                if (!IsFinite(transposeInverse) ||
                    std::abs(Determinant(transposeInverse)) <= MinPolarDeterminant)
                {
                    return Mat3::Identity();
                }

                const Mat3 nextRotation = Scale(Add(rotation, transposeInverse), 0.5f);
                if (!IsFinite(nextRotation))
                {
                    return Mat3::Identity();
                }

                if (FrobeniusDifferenceSquared(nextRotation, rotation) <= 1.0e-10f)
                {
                    rotation = nextRotation;
                    earlyExit = iterationCount < 12;
                    break;
                }

                rotation = nextRotation;
            }

            if (!IsFinite(rotation))
            {
                return Mat3::Identity();
            }

            if (Determinant(rotation) < 0.0f)
            {
                rotation.columns[2] *= -1.0f;
            }

            return rotation;
        }

        bool HasValidLocalPositions(const Tet& tet, const std::vector<Vec3>& positions)
        {
            const int nodeCount = static_cast<int>(positions.size());
            return tet.node0 >= 0 && tet.node0 < nodeCount &&
                tet.node1 >= 0 && tet.node1 < nodeCount &&
                tet.node2 >= 0 && tet.node2 < nodeCount &&
                tet.node3 >= 0 && tet.node3 < nodeCount;
        }

        Matrix6 BuildElasticityMatrix(float youngModulus, float poissonRatio)
        {
            Matrix6 d{};
            const float e = std::max(0.0f, youngModulus);
            if (e <= 0.0f)
            {
                return d;
            }

            const float nu = std::max(-0.99f, std::min(0.49f, poissonRatio));
            const float lambda = e * nu / ((1.0f + nu) * (1.0f - 2.0f * nu));
            const float mu = e / (2.0f * (1.0f + nu));

            d[0][0] = lambda + 2.0f * mu;
            d[0][1] = lambda;
            d[0][2] = lambda;
            d[1][0] = lambda;
            d[1][1] = lambda + 2.0f * mu;
            d[1][2] = lambda;
            d[2][0] = lambda;
            d[2][1] = lambda;
            d[2][2] = lambda + 2.0f * mu;
            d[3][3] = mu;
            d[4][4] = mu;
            d[5][5] = mu;
            return d;
        }

        Matrix6x12 BuildStrainDisplacementMatrix(const Tet& tet)
        {
            Matrix6x12 b{};

            for (int node = 0; node < 4; ++node)
            {
                const Vec3& grad = tet.shapeFunctionGradients[node];
                const int column = node * 3;

                // B maps nodal displacement to engineering small strain:
                // [exx, eyy, ezz, gamma_xy, gamma_yz, gamma_xz]^T = B * u_e.
                b[0][column + 0] = grad.x;
                b[1][column + 1] = grad.y;
                b[2][column + 2] = grad.z;
                b[3][column + 0] = grad.y;
                b[3][column + 1] = grad.x;
                b[4][column + 1] = grad.z;
                b[4][column + 2] = grad.y;
                b[5][column + 0] = grad.z;
                b[5][column + 2] = grad.x;
            }

            return b;
        }

        Vector12 BuildDisplacementVector(
            const Tet& tet,
            const std::vector<Vec3>& positions)
        {
            const int nodeIndices[4] = {tet.node0, tet.node1, tet.node2, tet.node3};
            Vector12 displacement{};

            for (int node = 0; node < 4; ++node)
            {
                const Vec3 u =
                    positions[static_cast<std::size_t>(nodeIndices[node])] -
                    tet.restPositions[node];
                displacement[node * 3 + 0] = u.x;
                displacement[node * 3 + 1] = u.y;
                displacement[node * 3 + 2] = u.z;
            }

            return displacement;
        }

        Vector12 BuildCorotatedDisplacementVector(
            const Tet& tet,
            const std::vector<Vec3>& positions,
            const Mat3& rotation)
        {
            const int nodeIndices[4] = {tet.node0, tet.node1, tet.node2, tet.node3};
            const Vec3& currentOrigin = positions[static_cast<std::size_t>(tet.node0)];
            const Vec3& restOrigin = tet.restPositions[0];
            const Mat3 inverseRotation = Transpose(rotation);
            Vector12 displacement{};

            for (int node = 0; node < 4; ++node)
            {
                const Vec3 localCurrent =
                    inverseRotation *
                    (positions[static_cast<std::size_t>(nodeIndices[node])] - currentOrigin) +
                    restOrigin;
                const Vec3 u = localCurrent - tet.restPositions[node];
                displacement[node * 3 + 0] = u.x;
                displacement[node * 3 + 1] = u.y;
                displacement[node * 3 + 2] = u.z;
            }

            return displacement;
        }

        Vector6 Multiply(const Matrix6x12& matrix, const Vector12& vector)
        {
            Vector6 result{};
            for (int row = 0; row < 6; ++row)
            {
                for (int column = 0; column < 12; ++column)
                {
                    result[row] += matrix[row][column] * vector[column];
                }
            }
            return result;
        }

        Vector6 Multiply(const Matrix6& matrix, const Vector6& vector)
        {
            Vector6 result{};
            for (int row = 0; row < 6; ++row)
            {
                for (int column = 0; column < 6; ++column)
                {
                    result[row] += matrix[row][column] * vector[column];
                }
            }
            return result;
        }

        Vector12 MultiplyTranspose(const Matrix6x12& matrix, const Vector6& vector)
        {
            Vector12 result{};
            for (int column = 0; column < 12; ++column)
            {
                for (int row = 0; row < 6; ++row)
                {
                    result[column] += matrix[row][column] * vector[row];
                }
            }
            return result;
        }

        Matrix12 BuildElementStiffness(const Matrix6x12& b, const Matrix6& d, float volume)
        {
            Matrix12 stiffness{};

            // K_e = V * B^T * D * B.
            for (int row = 0; row < 12; ++row)
            {
                for (int column = 0; column < 12; ++column)
                {
                    float value = 0.0f;
                    for (int i = 0; i < 6; ++i)
                    {
                        for (int j = 0; j < 6; ++j)
                        {
                            value += b[i][row] * d[i][j] * b[j][column];
                        }
                    }
                    stiffness[row][column] = volume * value;
                }
            }

            return stiffness;
        }

        void AssembleStiffness(
            const Tet& tet,
            const Matrix12& stiffness,
            TetElementContribution& contribution,
            const Mat3* rotation = nullptr)
        {
            const int nodeIndices[4] = { tet.node0, tet.node1, tet.node2, tet.node3 };
            const Mat3 rotationTranspose =
                rotation != nullptr ? Transpose(*rotation) : Mat3::Identity();

            for (int rowNode = 0; rowNode < 4; ++rowNode)
            {
                for (int columnNode = 0; columnNode < 4; ++columnNode)
                {
                    Mat3 block = Mat3::Zero();

                    for (int rowAxis = 0; rowAxis < 3; ++rowAxis)
                    {
                        for (int columnAxis = 0; columnAxis < 3; ++columnAxis)
                        {
                            SetBlockValue(
                                block,
                                rowAxis,
                                columnAxis,
                                stiffness[rowNode * 3 + rowAxis]
                                [columnNode * 3 + columnAxis]);
                        }
                    }

                    if (rotation != nullptr)
                    {
                        block = (*rotation) * block * rotationTranspose;
                    }

                    contribution.localNodeIndices[rowNode] = nodeIndices[rowNode];
                    contribution.localNodeIndices[columnNode] = nodeIndices[columnNode];
                    contribution.stiffness[rowNode][columnNode] = block;
                }
            }
        }

        TetElementContribution ComputeElasticForcesAndStiffness(
            const Tet& tet,
            const TetFemCache& cache,
            const std::vector<Vec3>& velocities,
            const Vector12& displacement,
            const Mat3* rotation)
        {
            TetElementContribution contribution;
            // epsilon = B * u_e.
            const Vector6 strain = Multiply(cache.B, displacement);

            // sigma = D * epsilon.
            const Vector6 stress = Multiply(cache.D, strain);

            // f_int = V * B^T * sigma. Return external force contribution -f_int.
            Vector12 internalForce = MultiplyTranspose(cache.B, stress);
            for (float& value : internalForce)
            {
                value *= tet.restVolume;
            }

            const int nodeIndices[4] = {tet.node0, tet.node1, tet.node2, tet.node3};
            for (int node = 0; node < 4; ++node)
            {
                const int nodeIndex = nodeIndices[node];
                contribution.localNodeIndices[node] = nodeIndex;
                Vec3 elasticForce{
                    -internalForce[node * 3 + 0],
                    -internalForce[node * 3 + 1],
                    -internalForce[node * 3 + 2]};
                if (rotation != nullptr)
                {
                    elasticForce = (*rotation) * elasticForce;
                }

                const Vec3 dampingForce =
                    velocities[static_cast<std::size_t>(nodeIndex)] * (-tet.damping);
                // Temporary per-node damping for stability.
                // This is not Rayleigh damping or element-level damping.
                // A future milestone should replace this with a proper damping model.
                contribution.forces[node] = elasticForce + dampingForce;
            }

            // Store positive element stiffness K_e.
            // The component-level assembly decides how to combine it into the solver.
            AssembleStiffness(tet, cache.Ke, contribution, rotation);
            return contribution;
        }
    }

    bool FEMModel::IsFemModelImplemented(FemModel femModel)
    {
        return femModel == FemModel::Linear ||
            femModel == FemModel::Corotational;
    }

    const char* FEMModel::GetNotImplementedMessage(FemModel femModel)
    {
        switch (femModel)
        {
        case FemModel::Corotational:
            return "";
        case FemModel::NeoHookean:
            return "NeoHookean FEM is not implemented yet";
        case FemModel::Linear:
            return "";
        }

        return "Unknown FEM model is not implemented";
    }

    bool FEMModel::ComputeForces(
        FemModel femModel,
        const std::vector<Tet>& tets,
        const std::vector<Vec3>& positions,
        const std::vector<Vec3>& velocities,
        std::vector<TetElementContribution>& outContributions)
    {
        std::vector<TetFemCache> tetFemCache;
        tetFemCache.reserve(tets.size());
        for (const Tet& tet : tets)
        {
            tetFemCache.push_back(BuildTetFemCache(tet));
        }

        return ComputeForces(
            femModel,
            tets,
            tetFemCache,
            positions,
            velocities,
            outContributions);
    }

    bool FEMModel::ComputeForces(
        FemModel femModel,
        const std::vector<Tet>& tets,
        const std::vector<TetFemCache>& tetFemCache,
        const std::vector<Vec3>& positions,
        const std::vector<Vec3>& velocities,
        std::vector<TetElementContribution>& outContributions)
    {
        switch (femModel)
        {
        case FemModel::Linear:
            ComputeElasticForces(
                tets,
                tetFemCache,
                positions,
                velocities,
                outContributions);
            return true;
        case FemModel::Corotational:
            ComputeCorotationalElasticForces(
                tets,
                tetFemCache,
                positions,
                velocities,
                outContributions);
            return true;
        case FemModel::NeoHookean:
            return false;
        }

        return false;
    }

    void FEMModel::InitializeTetRestData(
        Tet& tet,
        const std::vector<Vec3>& restPositions)
    {
        if (!HasValidLocalPositions(tet, restPositions))
        {
            return;
        }

        const Vec3& rest0 = restPositions[static_cast<std::size_t>(tet.node0)];
        const Vec3& rest1 = restPositions[static_cast<std::size_t>(tet.node1)];
        const Vec3& rest2 = restPositions[static_cast<std::size_t>(tet.node2)];
        const Vec3& rest3 = restPositions[static_cast<std::size_t>(tet.node3)];
        const Mat3 restDm = Mat3::FromColumns(rest1 - rest0, rest2 - rest0, rest3 - rest0);
        const float determinant = Determinant(restDm);
        tet.restPositions[0] = rest0;
        tet.restPositions[1] = rest1;
        tet.restPositions[2] = rest2;
        tet.restPositions[3] = rest3;

        const float volume = std::abs(determinant) / 6.0f;
        if (volume <= MinTetVolume)
        {
            tet.restVolume = 0.0f;
            tet.restDmInverse = Mat3::Zero();
            tet.shapeFunctionGradients[0] = Vec3{};
            tet.shapeFunctionGradients[1] = Vec3{};
            tet.shapeFunctionGradients[2] = Vec3{};
            tet.shapeFunctionGradients[3] = Vec3{};
            return;
        }

        tet.restVolume = volume;
        tet.restDmInverse = Inverse(restDm);

        // For linear tetrahedra, gradients of shape functions are constant in rest space.
        // DmInv maps world-space differential position to barycentric coordinates N1..N3.
        tet.shapeFunctionGradients[1] = GetRow(tet.restDmInverse, 0);
        tet.shapeFunctionGradients[2] = GetRow(tet.restDmInverse, 1);
        tet.shapeFunctionGradients[3] = GetRow(tet.restDmInverse, 2);
        tet.shapeFunctionGradients[0] =
            -(tet.shapeFunctionGradients[1] +
              tet.shapeFunctionGradients[2] +
              tet.shapeFunctionGradients[3]);
    }

    TetFemCache FEMModel::BuildTetFemCache(const Tet& tet)
    {
        TetFemCache cache;
        cache.B = BuildStrainDisplacementMatrix(tet);
        cache.D = BuildElasticityMatrix(tet.youngModulus, tet.poissonRatio);
        cache.Ke = BuildElementStiffness(cache.B, cache.D, tet.restVolume);
        return cache;
    }

    void FEMModel::ComputeElasticForces(
        const std::vector<Tet>& tets,
        const std::vector<Vec3>& positions,
        const std::vector<Vec3>& velocities,
        std::vector<TetElementContribution>& outContributions)
    {
        std::vector<TetFemCache> tetFemCache;
        tetFemCache.reserve(tets.size());
        for (const Tet& tet : tets)
        {
            tetFemCache.push_back(BuildTetFemCache(tet));
        }

        ComputeElasticForces(
            tets,
            tetFemCache,
            positions,
            velocities,
            outContributions);
    }

    void FEMModel::ComputeElasticForces(
        const std::vector<Tet>& tets,
        const std::vector<TetFemCache>& tetFemCache,
        const std::vector<Vec3>& positions,
        const std::vector<Vec3>& velocities,
        std::vector<TetElementContribution>& outContributions)
    {
        const std::size_t count = std::min(tets.size(), tetFemCache.size());
        for (std::size_t tetIndex = 0; tetIndex < count; ++tetIndex)
        {
            const Tet& tet = tets[tetIndex];
            if (!tet.active)
            {
                continue;
            }

            if (!HasValidLocalPositions(tet, positions) ||
                !HasValidLocalPositions(tet, velocities) ||
                tet.restVolume <= 0.0f)
            {
                continue;
            }

            const Vector12 displacement = BuildDisplacementVector(tet, positions);
            outContributions.push_back(
                ComputeElasticForcesAndStiffness(
                    tet,
                    tetFemCache[tetIndex],
                    velocities,
                    displacement,
                    nullptr));
        }
    }

    void FEMModel::ComputeCorotationalElasticForces(
        const std::vector<Tet>& tets,
        const std::vector<Vec3>& positions,
        const std::vector<Vec3>& velocities,
        std::vector<TetElementContribution>& outContributions)
    {
        std::vector<TetFemCache> tetFemCache;
        tetFemCache.reserve(tets.size());
        for (const Tet& tet : tets)
        {
            tetFemCache.push_back(BuildTetFemCache(tet));
        }

        ComputeCorotationalElasticForces(
            tets,
            tetFemCache,
            positions,
            velocities,
            outContributions);
    }

    void FEMModel::ComputeCorotationalElasticForces(
        const std::vector<Tet>& tets,
        const std::vector<TetFemCache>& tetFemCache,
        const std::vector<Vec3>& positions,
        const std::vector<Vec3>& velocities,
        std::vector<TetElementContribution>& outContributions)
    {
        const std::size_t count = std::min(tets.size(), tetFemCache.size());
        for (std::size_t tetIndex = 0; tetIndex < count; ++tetIndex)
        {
            const Tet& tet = tets[tetIndex];
            if (!tet.active)
            {
                continue;
            }

            if (!HasValidLocalPositions(tet, positions) ||
                !HasValidLocalPositions(tet, velocities) ||
                tet.restVolume <= 0.0f)
            {
                continue;
            }

            Mat3 deformationGradient;
            {
                const Mat3 ds = BuildDm(tet, positions);
                deformationGradient = ds * tet.restDmInverse;
            }
            if (!IsFinite(deformationGradient))
            {
                continue;
            }

            const Mat3 rotation = ExtractRotationPolar(deformationGradient);
            const Vector12 displacement =
                BuildCorotatedDisplacementVector(tet, positions, rotation);
            outContributions.push_back(
                ComputeElasticForcesAndStiffness(
                    tet,
                    tetFemCache[tetIndex],
                    velocities,
                    displacement,
                    &rotation));
        }
    }

    void FEMModel::ComputeLumpedMass(
        const Material& material,
        const std::vector<Tet>& tets,
        std::vector<TetMassContribution>& outContributions)
    {
        const float density = std::max(0.0f, material.density);
        if (!std::isfinite(density))
        {
            return;
        }

        for (const Tet& tet : tets)
        {
            if (!tet.active)
            {
                continue;
            }

            if (!std::isfinite(tet.restVolume) || tet.restVolume <= 0.0f)
            {
                continue;
            }

            const float nodalMass = density * tet.restVolume * 0.25f;
            TetMassContribution contribution;
            contribution.localNodeIndices[0] = tet.node0;
            contribution.localNodeIndices[1] = tet.node1;
            contribution.localNodeIndices[2] = tet.node2;
            contribution.localNodeIndices[3] = tet.node3;
            contribution.nodalMass = nodalMass;
            outContributions.push_back(contribution);
        }
    }

    std::vector<std::pair<int, int>> FEMModel::BuildSparsePatternFromTetConnectivity(
        const std::vector<Tet>& tets)
    {
        std::vector<std::pair<int, int>> blockCoordinates;
        blockCoordinates.reserve(tets.size() * 16u);

        for (const Tet& tet : tets)
        {
            const int nodes[4] = { tet.node0, tet.node1, tet.node2, tet.node3 };
            for (int row = 0; row < 4; ++row)
            {
                for (int column = 0; column < 4; ++column)
                {
                    blockCoordinates.push_back({ nodes[row], nodes[column] });
                }
            }
        }

        return blockCoordinates;
    }
}
