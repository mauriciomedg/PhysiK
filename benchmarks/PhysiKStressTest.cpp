#include "PhysiK/Core/Solvers/Linear/LinearSolver.h"
#include "PhysiK/Core/Solvers/SolverData.h"
#include "PhysiK/Math/SparseBlockMatrix.h"
#include "PhysiK/Math/Vec3.h"
#include "PhysiK/PhysicsData/Node.h"

#include <chrono>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <vector>

namespace
{
    struct StressSystem
    {
        std::vector<PhysiK::Node> nodes;
        PhysiK::SolverData solverData;
    };

    struct TimingStats
    {
        double total = 0.0;
        double minimum = std::numeric_limits<double>::max();
        double maximum = 0.0;
        int count = 0;

        void Add(double value)
        {
            total += value;
            minimum = std::min(minimum, value);
            maximum = std::max(maximum, value);
            ++count;
        }

        double Average() const
        {
            return count > 0 ? total / static_cast<double>(count) : 0.0;
        }

        double Minimum() const
        {
            return count > 0 ? minimum : 0.0;
        }
    };

#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
    struct ProfileStats
    {
        TimingStats totalSolve;
        TimingStats sparseMultiply;
        TimingStats dotProduct;
        TimingStats vectorUpdate;
        TimingStats preconditionerSetup;
        TimingStats preconditionerApply;
        TimingStats iterations;
        TimingStats residualNorm;

        void Add(const PhysiK::LinearSolverProfileData& profile)
        {
            totalSolve.Add(profile.totalSolveMilliseconds);
            sparseMultiply.Add(profile.sparseMatrixMultiplyMilliseconds);
            dotProduct.Add(profile.dotProductMilliseconds);
            vectorUpdate.Add(profile.vectorUpdateMilliseconds);
            preconditionerSetup.Add(profile.preconditionerSetupMilliseconds);
            preconditionerApply.Add(profile.preconditionerApplyMilliseconds);
            iterations.Add(static_cast<double>(profile.iterations));
            residualNorm.Add(static_cast<double>(profile.residualNorm));
        }
    };
#endif

    void Fail(const char* message)
    {
        std::fprintf(stderr, "PhysiKStressTest failed: %s\n", message);
        std::exit(1);
    }

    PhysiK::Mat3 DiagonalBlock(float value)
    {
        return PhysiK::Mat3::FromColumns(
            PhysiK::Vec3{value, 0.0f, 0.0f},
            PhysiK::Vec3{0.0f, value, 0.0f},
            PhysiK::Vec3{0.0f, 0.0f, value});
    }

    int NodeIndex(int x, int y, int z, int width, int height)
    {
        return x + width * (y + height * z);
    }

    StressSystem BuildStructuredImplicitFemLikeSystem(int width, int height, int depth)
    {
        StressSystem system;
        const int nodeCount = width * height * depth;
        system.nodes.resize(static_cast<std::size_t>(nodeCount));

        for (int z = 0; z < depth; ++z)
        {
            for (int y = 0; y < height; ++y)
            {
                for (int x = 0; x < width; ++x)
                {
                    const int nodeIndex = NodeIndex(x, y, z, width, height);
                    PhysiK::Node& node = system.nodes[static_cast<std::size_t>(nodeIndex)];
                    node.position = PhysiK::Vec3{
                        static_cast<float>(x),
                        static_cast<float>(y),
                        static_cast<float>(z)};
                    node.fixed = x == 0;

                    if (!node.fixed)
                    {
                        system.solverData.AddNodeMass(nodeIndex, 1.0f);
                        system.solverData.AddNodeForce(
                            nodeIndex,
                            PhysiK::Vec3{
                                0.02f * static_cast<float>((nodeIndex % 7) + 1),
                                -0.03f,
                                0.01f * static_cast<float>((nodeIndex % 5) + 1)});
                    }
                }
            }
        }

        auto addEdge = [&](int a, int b)
        {
            constexpr float stiffness = 120.0f;
            system.solverData.AddStiffnessBlock(a, a, DiagonalBlock(stiffness));
            system.solverData.AddStiffnessBlock(a, b, DiagonalBlock(-stiffness));
            system.solverData.AddStiffnessBlock(b, a, DiagonalBlock(-stiffness));
            system.solverData.AddStiffnessBlock(b, b, DiagonalBlock(stiffness));
        };

        for (int z = 0; z < depth; ++z)
        {
            for (int y = 0; y < height; ++y)
            {
                for (int x = 0; x < width; ++x)
                {
                    const int node = NodeIndex(x, y, z, width, height);
                    if (x + 1 < width)
                    {
                        addEdge(node, NodeIndex(x + 1, y, z, width, height));
                    }
                    if (y + 1 < height)
                    {
                        addEdge(node, NodeIndex(x, y + 1, z, width, height));
                    }
                    if (z + 1 < depth)
                    {
                        addEdge(node, NodeIndex(x, y, z + 1, width, height));
                    }
                }
            }
        }

        return system;
    }

    double DeltaVelocityNorm(const std::vector<float>& values)
    {
        double sum = 0.0;
        for (float value : values)
        {
            sum += static_cast<double>(value) * static_cast<double>(value);
        }
        return std::sqrt(sum);
    }

    bool SolveOnce(
        const std::vector<PhysiK::Node>& nodes,
        PhysiK::SolverData solverData,
        double& outMilliseconds,
        double& outDeltaVelocityNorm)
    {
        constexpr float dt = 0.016f;
        if (!solverData.PrecomputeImplicitSolve(nodes, dt))
        {
            Fail("could not precompute implicit solve");
        }

        const auto start = std::chrono::steady_clock::now();
        const bool solved = solverData.SolveImplicitLinearSystem();
        const auto end = std::chrono::steady_clock::now();
        if (!solved)
        {
            return false;
        }

        outMilliseconds =
            std::chrono::duration<double, std::milli>(end - start).count();
        outDeltaVelocityNorm = DeltaVelocityNorm(solverData.GetDeltaVelocity());
        return true;
    }

    int ParsePositiveInt(const char* value, int fallback)
    {
        if (value == nullptr)
        {
            return fallback;
        }

        const int parsed = std::atoi(value);
        return parsed > 0 ? parsed : fallback;
    }

    PhysiK::SparseBlockMatrixMultiplyMode ParseMultiplyMode(const char* value)
    {
        if (value == nullptr || std::strcmp(value, "serial") == 0)
        {
            return PhysiK::SparseBlockMatrixMultiplyMode::Serial;
        }

        if (std::strcmp(value, "parallel") == 0 ||
            std::strcmp(value, "cv") == 0)
        {
            return PhysiK::SparseBlockMatrixMultiplyMode::ConditionVariableParallel;
        }

        if (std::strcmp(value, "spin") == 0 ||
            std::strcmp(value, "spinning") == 0)
        {
            return PhysiK::SparseBlockMatrixMultiplyMode::SpinningWorkers;
        }

        return PhysiK::SparseBlockMatrixMultiplyMode::Serial;
    }

    const char* GetMultiplyModeName(PhysiK::SparseBlockMatrixMultiplyMode mode)
    {
        switch (mode)
        {
        case PhysiK::SparseBlockMatrixMultiplyMode::ConditionVariableParallel:
            return "parallel";
        case PhysiK::SparseBlockMatrixMultiplyMode::SpinningWorkers:
            return "spin";
        case PhysiK::SparseBlockMatrixMultiplyMode::Serial:
        default:
            return "serial";
        }
    }

    void PrintStats(const char* label, const TimingStats& stats, const char* unit)
    {
        if (unit[0] == '\0')
        {
            std::printf(
                "  %-28s avg=%12.8g  min=%12.8g  max=%12.8g\n",
                label,
                stats.Average(),
                stats.Minimum(),
                stats.maximum);
            return;
        }

        std::printf(
            "  %-28s avg=%9.4f %s  min=%9.4f %s  max=%9.4f %s\n",
            label,
            stats.Average(),
            unit,
            stats.Minimum(),
            unit,
            stats.maximum,
            unit);
    }
}

int main(int argc, char** argv)
{
    const int width = ParsePositiveInt(argc > 1 ? argv[1] : nullptr, 20);
    const int height = ParsePositiveInt(argc > 2 ? argv[2] : nullptr, 20);
    const int depth = ParsePositiveInt(argc > 3 ? argv[3] : nullptr, 10);
    const int iterations = ParsePositiveInt(argc > 4 ? argv[4] : nullptr, 200);
    const PhysiK::SparseBlockMatrixMultiplyMode multiplyMode =
        ParseMultiplyMode(argc > 5 ? argv[5] : nullptr);
    const int warmups = std::min(10, std::max(1, iterations / 10));
    PhysiK::SetSparseBlockMatrixMultiplyMode(multiplyMode);

    std::printf(
        "PhysiK stress test grid: %d x %d x %d nodes\n",
        width,
        height,
        depth);
    std::printf(
        "Warmups: %d, benchmark iterations: %d\n",
        warmups,
        iterations);
    std::printf(
        "SparseBlockMatrix multiply mode: %s\n",
        GetMultiplyModeName(multiplyMode));

    const StressSystem system = BuildStructuredImplicitFemLikeSystem(width, height, depth);
    for (int i = 0; i < warmups; ++i)
    {
        double milliseconds = 0.0;
        double norm = 0.0;
        if (!SolveOnce(system.nodes, system.solverData, milliseconds, norm))
        {
            std::fprintf(stderr, "CurrentLinearSolver warmup did not converge\n");
            return 1;
        }
    }

    TimingStats wallClockStats;
    TimingStats deltaVelocityNormStats;
#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
    ProfileStats profileStats;
#endif

    for (int i = 0; i < iterations; ++i)
    {
        double milliseconds = 0.0;
        double norm = 0.0;
        if (!SolveOnce(system.nodes, system.solverData, milliseconds, norm))
        {
            std::fprintf(stderr, "CurrentLinearSolver benchmark solve did not converge\n");
            return 1;
        }

        wallClockStats.Add(milliseconds);
        deltaVelocityNormStats.Add(norm);
#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
        profileStats.Add(PhysiK::GetCurrentLinearSolverProfile());
#endif
    }

    std::printf("CurrentLinearSolver wall-clock solve timing:\n");
    std::printf(
        "  SparseBlockMatrix worker count: %d\n",
        PhysiK::GetSparseBlockMatrixLastMultiplyThreadCount());
    PrintStats("solve", wallClockStats, "ms");
    PrintStats("delta velocity norm", deltaVelocityNormStats, "");

#if defined(PHYSIK_ENABLE_SOLVER_PROFILING)
    std::printf("CurrentLinearSolver internal profiling:\n");
    PrintStats("total solve", profileStats.totalSolve, "ms");
    PrintStats("SparseBlockMatrix multiply", profileStats.sparseMultiply, "ms");
    PrintStats("dot product", profileStats.dotProduct, "ms");
    PrintStats("vector update", profileStats.vectorUpdate, "ms");
    PrintStats("preconditioner setup", profileStats.preconditionerSetup, "ms");
    PrintStats("preconditioner apply", profileStats.preconditionerApply, "ms");
    PrintStats("iterations", profileStats.iterations, "");
    PrintStats("residual norm", profileStats.residualNorm, "");
#else
    std::printf("Internal profiling is disabled in this build.\n");
#endif

    return 0;
}
