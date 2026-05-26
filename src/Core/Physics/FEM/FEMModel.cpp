#include "PhysiK/Core/Physics/FEM/FEMModel.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <optional>

#include "PhysiK/Components/TetMeshComponent.h"
#if defined(PHYSIK_ENABLE_PERF_LOGGING)
#include "PhysiK/Core/Performance/PerformanceLogger.h"
#endif
#include "PhysiK/Core/Solvers/SolverData.h"
#include "PhysiK/Core/World/World.h"

namespace PhysiK
{
    namespace
    {
        using Vector6 = std::array<float, 6>;
        using Vector12 = std::array<float, 12>;

        constexpr float MinTetVolume = 1.0e-8f;
        constexpr float MinPolarDeterminant = 1.0e-8f;

#if defined(PHYSIK_ENABLE_PERF_LOGGING)
        thread_local PerformanceLogRecord* currentPerformanceRecord = nullptr;


        void AddElapsed(double& target, const PerformanceTimer& timer)
        {
            target += timer.ElapsedMilliseconds();
        }

        void RecordPolarCall(int iterations, bool earlyExit)
        {
            if (currentPerformanceRecord == nullptr)
            {
                return;
            }

            ++currentPerformanceRecord->polarCallCount;
            const int callCount = currentPerformanceRecord->polarCallCount;
            currentPerformanceRecord->averagePolarIterations +=
                (static_cast<double>(iterations) -
                 currentPerformanceRecord->averagePolarIterations) /
                static_cast<double>(callCount);
            currentPerformanceRecord->maxPolarIterationsObserved = std::max(
                currentPerformanceRecord->maxPolarIterationsObserved,
                iterations);
            if (earlyExit)
            {
                ++currentPerformanceRecord->polarEarlyExitCount;
            }
        }

        void AddPolarElapsed(double elapsedMilliseconds)
        {
            if (currentPerformanceRecord == nullptr)
            {
                return;
            }

            currentPerformanceRecord->extractRotationPolarMs += elapsedMilliseconds;
            const int callCount = currentPerformanceRecord->polarCallCount;
            if (callCount > 0)
            {
                currentPerformanceRecord->averageExtractRotationPolarMs +=
                    (elapsedMilliseconds -
                     currentPerformanceRecord->averageExtractRotationPolarMs) /
                    static_cast<double>(callCount);
            }
        }
#endif

        Mat3 BuildDm(const Tet& tet, const std::vector<Node>& nodes)
        {
            const Vec3& x0 = nodes[static_cast<std::size_t>(tet.node0)].position;
            const Vec3& x1 = nodes[static_cast<std::size_t>(tet.node1)].position;
            const Vec3& x2 = nodes[static_cast<std::size_t>(tet.node2)].position;
            const Vec3& x3 = nodes[static_cast<std::size_t>(tet.node3)].position;
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
#if defined(PHYSIK_ENABLE_PERF_LOGGING)
                RecordPolarCall(0, true);
#endif
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
#if defined(PHYSIK_ENABLE_PERF_LOGGING)
                    RecordPolarCall(iterationCount, true);
#endif
                    return Mat3::Identity();
                }

                const Mat3 nextRotation = Scale(Add(rotation, transposeInverse), 0.5f);
                if (!IsFinite(nextRotation))
                {
#if defined(PHYSIK_ENABLE_PERF_LOGGING)
                    RecordPolarCall(iterationCount, true);
#endif
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
#if defined(PHYSIK_ENABLE_PERF_LOGGING)
                RecordPolarCall(iterationCount, true);
#endif
                return Mat3::Identity();
            }

            if (Determinant(rotation) < 0.0f)
            {
                rotation.columns[2] *= -1.0f;
            }

#if defined(PHYSIK_ENABLE_PERF_LOGGING)
            RecordPolarCall(iterationCount, earlyExit);
#endif
            return rotation;
        }

        bool HasValidNodes(const Tet& tet, const std::vector<Node>& nodes)
        {
            const int nodeCount = static_cast<int>(nodes.size());
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

        Vector12 BuildDisplacementVector(const Tet& tet, const std::vector<Node>& nodes)
        {
            const int nodeIndices[4] = {tet.node0, tet.node1, tet.node2, tet.node3};
            Vector12 displacement{};

            for (int node = 0; node < 4; ++node)
            {
                const Vec3 u =
                    nodes[static_cast<std::size_t>(nodeIndices[node])].position -
                    tet.restPositions[node];
                displacement[node * 3 + 0] = u.x;
                displacement[node * 3 + 1] = u.y;
                displacement[node * 3 + 2] = u.z;
            }

            return displacement;
        }

        Vector12 BuildCorotatedDisplacementVector(
            const Tet& tet,
            const std::vector<Node>& nodes,
            const Mat3& rotation)
        {
            const int nodeIndices[4] = {tet.node0, tet.node1, tet.node2, tet.node3};
            const Vec3& currentOrigin = nodes[static_cast<std::size_t>(tet.node0)].position;
            const Vec3& restOrigin = tet.restPositions[0];
            const Mat3 inverseRotation = Transpose(rotation);
            Vector12 displacement{};

            for (int node = 0; node < 4; ++node)
            {
                const Vec3 localCurrent =
                    inverseRotation *
                    (nodes[static_cast<std::size_t>(nodeIndices[node])].position - currentOrigin) +
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
            SolverData& solverData,
            const Mat3* rotation = nullptr)
        {
            const int nodeIndices[4] = { tet.node0, tet.node1, tet.node2, tet.node3 };
            const Mat3 rotationTranspose =
                rotation != nullptr ? Transpose(*rotation) : Mat3::Identity();

#if defined(PHYSIK_ENABLE_PERF_LOGGING)
            PerformanceLogRecord* performanceRecord = currentPerformanceRecord;
            const bool logPerformance = performanceRecord != nullptr;

            std::optional<PerformanceTimer> assembleStiffnessTimer;
            if (logPerformance)
            {
                assembleStiffnessTimer.emplace();
            }
#endif

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

                    solverData.AddStiffnessBlock(
                        nodeIndices[rowNode],
                        nodeIndices[columnNode],
                        block);
                }
            }

#if defined(PHYSIK_ENABLE_PERF_LOGGING)
            if (logPerformance)
            {
                AddElapsed(
                    performanceRecord->assembleStiffnessBlocksMs,
                    *assembleStiffnessTimer);
            }
#endif
        }

        void AddElasticForcesAndStiffness(
            const Tet& tet,
            const TetFemCache& cache,
            const std::vector<Node>& nodes,
            SolverData& solverData,
            const Vector12& displacement,
            const Mat3* rotation)
        {
#if defined(PHYSIK_ENABLE_PERF_LOGGING)
            PerformanceLogRecord* performanceRecord = currentPerformanceRecord;
            const bool logPerformance = performanceRecord != nullptr;
            std::optional<PerformanceTimer> forceTimer;
            if (logPerformance)
            {
                forceTimer.emplace();
            }
#endif
            // epsilon = B * u_e.
            const Vector6 strain = Multiply(cache.B, displacement);

            // sigma = D * epsilon.
            const Vector6 stress = Multiply(cache.D, strain);

            // f_int = V * B^T * sigma. SolverData stores total force, so assemble -f_int.
            Vector12 internalForce = MultiplyTranspose(cache.B, stress);
            for (float& value : internalForce)
            {
                value *= tet.restVolume;
            }
#if defined(PHYSIK_ENABLE_PERF_LOGGING)
            if (logPerformance)
            {
                const double forceMilliseconds = forceTimer->ElapsedMilliseconds();
                performanceRecord->assembleTetForceMs += forceMilliseconds;
                performanceRecord->computeElasticForcesMs += forceMilliseconds;
            }

            std::optional<PerformanceTimer> rhsWriteTimer;
            if (logPerformance)
            {
                rhsWriteTimer.emplace();
            }
#endif

            const int nodeIndices[4] = {tet.node0, tet.node1, tet.node2, tet.node3};
            for (int node = 0; node < 4; ++node)
            {
                const int nodeIndex = nodeIndices[node];
                Vec3 elasticForce{
                    -internalForce[node * 3 + 0],
                    -internalForce[node * 3 + 1],
                    -internalForce[node * 3 + 2]};
                if (rotation != nullptr)
                {
                    elasticForce = (*rotation) * elasticForce;
                }

                const Vec3 dampingForce =
                    nodes[static_cast<std::size_t>(nodeIndex)].velocity * (-tet.damping);
                // Temporary per-node damping for stability.
                // This is not Rayleigh damping or element-level damping.
                // A future milestone should replace this with a proper damping model.
                solverData.AddNodeForce(nodeIndex, elasticForce + dampingForce);
            }
#if defined(PHYSIK_ENABLE_PERF_LOGGING)
            if (logPerformance)
            {
                AddElapsed(performanceRecord->assembleRhsWriteMs, *rhsWriteTimer);
            }

            std::optional<PerformanceTimer> stiffnessTimer;
            if (logPerformance)
            {
                stiffnessTimer.emplace();
            }
#endif

            // Store positive element stiffness K_e.
            // The implicit solver decides how to combine it into the global system.
            AssembleStiffness(tet, cache.Ke, solverData, rotation);
#if defined(PHYSIK_ENABLE_PERF_LOGGING)
            if (logPerformance)
            {
                AddElapsed(performanceRecord->assembleTetStiffnessMs, *stiffnessTimer);
            }
#endif
        }
    }

    void FEMModel::UpdateSystem(
        World& world,
        TetMeshPhysicsComponent& owner,
        SolverData& solverData,
        float dt)
    {
        (void)dt;
        AccumulateForces(
            owner.GetFemModel(),
            owner.worldTets,
            owner.tetFemCache,
            world.GetNodes(),
            solverData);
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

    bool FEMModel::AccumulateForces(
        FemModel femModel,
        const std::vector<Tet>& tets,
        const std::vector<Node>& nodes,
        SolverData& solverData)
    {
        std::vector<TetFemCache> tetFemCache;
        tetFemCache.reserve(tets.size());
        for (const Tet& tet : tets)
        {
            tetFemCache.push_back(BuildTetFemCache(tet));
        }

        return AccumulateForces(femModel, tets, tetFemCache, nodes, solverData);
    }

    bool FEMModel::AccumulateForces(
        FemModel femModel,
        const std::vector<Tet>& tets,
        const std::vector<TetFemCache>& tetFemCache,
        const std::vector<Node>& nodes,
        SolverData& solverData)
    {
        switch (femModel)
        {
        case FemModel::Linear:
            AccumulateElasticForces(tets, tetFemCache, nodes, solverData);
            return true;
        case FemModel::Corotational:
            AccumulateCorotationalElasticForces(tets, tetFemCache, nodes, solverData);
            return true;
        case FemModel::NeoHookean:
            std::cerr << GetNotImplementedMessage(femModel) << '\n';
            return false;
        }

        std::cerr << GetNotImplementedMessage(femModel) << '\n';
        return false;
    }

#if defined(PHYSIK_ENABLE_PERF_LOGGING)
    void FEMModel::SetPerformanceLogRecord(PerformanceLogRecord* record)
    {
        currentPerformanceRecord = record;
    }

    PerformanceLogRecord* FEMModel::GetPerformanceLogRecord()
    {
        return currentPerformanceRecord;
    }
#endif

    void FEMModel::InitializeTetRestData(Tet& tet, const std::vector<Node>& nodes)
    {
        if (!HasValidNodes(tet, nodes))
        {
            return;
        }

        const Mat3 restDm = BuildDm(tet, nodes);
        const float determinant = Determinant(restDm);
        tet.restPositions[0] = nodes[static_cast<std::size_t>(tet.node0)].position;
        tet.restPositions[1] = nodes[static_cast<std::size_t>(tet.node1)].position;
        tet.restPositions[2] = nodes[static_cast<std::size_t>(tet.node2)].position;
        tet.restPositions[3] = nodes[static_cast<std::size_t>(tet.node3)].position;

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

    void FEMModel::AccumulateElasticForces(
        const std::vector<Tet>& tets,
        const std::vector<Node>& nodes,
        SolverData& solverData)
    {
        std::vector<TetFemCache> tetFemCache;
        tetFemCache.reserve(tets.size());
        for (const Tet& tet : tets)
        {
            tetFemCache.push_back(BuildTetFemCache(tet));
        }

        AccumulateElasticForces(tets, tetFemCache, nodes, solverData);
    }

    void FEMModel::AccumulateElasticForces(
        const std::vector<Tet>& tets,
        const std::vector<TetFemCache>& tetFemCache,
        const std::vector<Node>& nodes,
        SolverData& solverData)
    {
#if defined(PHYSIK_ENABLE_PERF_LOGGING)
        PerformanceLogRecord* performanceRecord = currentPerformanceRecord;
        const bool logPerformance = performanceRecord != nullptr;
        std::optional<PerformanceTimer> linearTimer;
        if (logPerformance)
        {
            linearTimer.emplace();
        }
#endif
        const std::size_t count = std::min(tets.size(), tetFemCache.size());
        for (std::size_t tetIndex = 0; tetIndex < count; ++tetIndex)
        {
            const Tet& tet = tets[tetIndex];
            if (!tet.active)
            {
                continue;
            }

            if (!HasValidNodes(tet, nodes) || tet.restVolume <= 0.0f)
            {
                continue;
            }

#if defined(PHYSIK_ENABLE_PERF_LOGGING)
            if (logPerformance)
            {
                ++performanceRecord->assembleTetCount;
            }
#endif
            const Vector12 displacement = BuildDisplacementVector(tet, nodes);
            AddElasticForcesAndStiffness(
                tet,
                tetFemCache[tetIndex],
                nodes,
                solverData,
                displacement,
                nullptr);
        }
#if defined(PHYSIK_ENABLE_PERF_LOGGING)
        if (logPerformance)
        {
            AddElapsed(performanceRecord->assembleLinearFemMs, *linearTimer);
        }
#endif
    }

    void FEMModel::AccumulateCorotationalElasticForces(
        const std::vector<Tet>& tets,
        const std::vector<Node>& nodes,
        SolverData& solverData)
    {
        std::vector<TetFemCache> tetFemCache;
        tetFemCache.reserve(tets.size());
        for (const Tet& tet : tets)
        {
            tetFemCache.push_back(BuildTetFemCache(tet));
        }

        AccumulateCorotationalElasticForces(tets, tetFemCache, nodes, solverData);
    }

    void FEMModel::AccumulateCorotationalElasticForces(
        const std::vector<Tet>& tets,
        const std::vector<TetFemCache>& tetFemCache,
        const std::vector<Node>& nodes,
        SolverData& solverData)
    {
#if defined(PHYSIK_ENABLE_PERF_LOGGING)
        PerformanceLogRecord* performanceRecord = currentPerformanceRecord;
        const bool logPerformance = performanceRecord != nullptr;
        std::optional<PerformanceTimer> corotationalTimer;
        if (logPerformance)
        {
            corotationalTimer.emplace();
        }
#endif
        const std::size_t count = std::min(tets.size(), tetFemCache.size());
        for (std::size_t tetIndex = 0; tetIndex < count; ++tetIndex)
        {
            const Tet& tet = tets[tetIndex];
            if (!tet.active)
            {
                continue;
            }

            if (!HasValidNodes(tet, nodes) || tet.restVolume <= 0.0f)
            {
                continue;
            }

#if defined(PHYSIK_ENABLE_PERF_LOGGING)
            if (logPerformance)
            {
                ++performanceRecord->assembleTetCount;
            }
#endif
            Mat3 deformationGradient;
#if defined(PHYSIK_ENABLE_PERF_LOGGING)
            std::optional<PerformanceTimer> deformationGradientTimer;
            if (logPerformance)
            {
                deformationGradientTimer.emplace();
            }
#endif
            {
                const Mat3 ds = BuildDm(tet, nodes);
                deformationGradient = ds * tet.restDmInverse;
            }
#if defined(PHYSIK_ENABLE_PERF_LOGGING)
            if (logPerformance)
            {
                AddElapsed(
                    performanceRecord->computeDeformationGradientMs,
                    *deformationGradientTimer);
            }
#endif
            if (!IsFinite(deformationGradient))
            {
                continue;
            }

#if defined(PHYSIK_ENABLE_PERF_LOGGING)
            std::optional<PerformanceTimer> polarTimer;
            if (logPerformance)
            {
                polarTimer.emplace();
            }
#endif
            const Mat3 rotation = ExtractRotationPolar(deformationGradient);
#if defined(PHYSIK_ENABLE_PERF_LOGGING)
            if (logPerformance)
            {
                AddPolarElapsed(polarTimer->ElapsedMilliseconds());
            }
#endif
            const Vector12 displacement = BuildCorotatedDisplacementVector(tet, nodes, rotation);
            AddElasticForcesAndStiffness(
                tet,
                tetFemCache[tetIndex],
                nodes,
                solverData,
                displacement,
                &rotation);
        }
#if defined(PHYSIK_ENABLE_PERF_LOGGING)
        if (logPerformance)
        {
            AddElapsed(performanceRecord->assembleCorotationalFemMs, *corotationalTimer);
        }
#endif
    }

    void FEMModel::AddLumpedMassToSolverData(
        const World& world,
        SolverData& solverData,
        int nodeIndex,
        float mass)
    {
        if (nodeIndex < 0 || !std::isfinite(mass))
        {
            return;
        }

        if (world.IsNodeFixed(nodeIndex))
        {
            return;
        }

        solverData.AddNodeMass(nodeIndex, std::max(0.0f, mass));
    }

    void FEMModel::AssembleLumpedMass(
        const TetMeshPhysicsComponent& component,
        const World& world,
        SolverData& solverData)
    {
        const float density = std::max(0.0f, component.material.density);
        if (!std::isfinite(density))
        {
            return;
        }

        for (const Tet& tet : component.worldTets)
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
            AddLumpedMassToSolverData(world, solverData, tet.node0, nodalMass);
            AddLumpedMassToSolverData(world, solverData, tet.node1, nodalMass);
            AddLumpedMassToSolverData(world, solverData, tet.node2, nodalMass);
            AddLumpedMassToSolverData(world, solverData, tet.node3, nodalMass);
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
