#include "PhysiK/Core/Solvers/Linear/LinearSolver.h"
#include "PhysiK/Core/Solvers/SolverData.h"
#include "PhysiK/Math/Vec3.h"
#include "PhysiK/PhysicsData/Node.h"

#include <chrono>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
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
                    node.stateIndex = nodeIndex;
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

    double DeltaVelocityNorm(const std::vector<PhysiK::Vec3>& values)
    {
        double sum = 0.0;
        for (const PhysiK::Vec3& value : values)
        {
            sum += static_cast<double>(value.x) * static_cast<double>(value.x);
            sum += static_cast<double>(value.y) * static_cast<double>(value.y);
            sum += static_cast<double>(value.z) * static_cast<double>(value.z);
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
        std::vector<PhysiK::Vec3> nodeVelocities(nodes.size(), PhysiK::Vec3{});
        if (!solverData.PrecomputeImplicitSolve(nodes, nodeVelocities, dt))
        {
            Fail("could not precompute implicit solve");
        }

        const auto start = std::chrono::steady_clock::now();
        PhysiK::ConjugateGradientSettings settings;
        const bool solved = solverData.SolveImplicitLinearSystem(settings);
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
    const int warmups = std::min(10, std::max(1, iterations / 10));

    std::printf(
        "PhysiK stress test grid: %d x %d x %d nodes\n",
        width,
        height,
        depth);
    std::printf(
        "Warmups: %d, benchmark iterations: %d\n",
        warmups,
        iterations);

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
    }

    std::printf("CurrentLinearSolver wall-clock solve timing:\n");
    PrintStats("solve", wallClockStats, "ms");
    PrintStats("delta velocity norm", deltaVelocityNormStats, "");

    return 0;
}
