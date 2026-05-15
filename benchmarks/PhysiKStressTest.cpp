#include "PhysiK/Core/Solvers/Linear/LinearSolver.h"
#include "PhysiK/Core/Solvers/SolverData.h"
#include "PhysiK/Math/Vec3.h"
#include "PhysiK/PhysicsData/Node.h"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace
{
    struct StressSystem
    {
        std::vector<PhysiK::Node> nodes;
        PhysiK::SolverData solverData;
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

    double SolveAndTime(
        const char* label,
        PhysiK::LinearSolverBackend backend,
        const std::vector<PhysiK::Node>& nodes,
        PhysiK::SolverData solverData)
    {
        constexpr float dt = 0.016f;
        if (!solverData.PrecomputeImplicitSolve(nodes, dt))
        {
            Fail("could not precompute implicit solve");
        }

        solverData.SetLinearSolverBackend(backend);
        const auto start = std::chrono::steady_clock::now();
        const bool solved = solverData.SolveImplicitLinearSystem();
        const auto end = std::chrono::steady_clock::now();
        if (!solved)
        {
            std::fprintf(stderr, "%s backend did not converge\n", label);
            return -1.0;
        }

        const double milliseconds =
            std::chrono::duration<double, std::milli>(end - start).count();
        const double norm = DeltaVelocityNorm(solverData.GetDeltaVelocity());
        std::printf(
            "%s: %.3f ms, dynamic blocks=%d, delta velocity norm=%.6f\n",
            label,
            milliseconds,
            solverData.GetDynamicBlockCount(),
            norm);
        return milliseconds;
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
}

int main(int argc, char** argv)
{
    const int width = ParsePositiveInt(argc > 1 ? argv[1] : nullptr, 6);
    const int height = ParsePositiveInt(argc > 2 ? argv[2] : nullptr, 6);
    const int depth = ParsePositiveInt(argc > 3 ? argv[3] : nullptr, 3);

    std::printf(
        "PhysiK stress test grid: %d x %d x %d nodes\n",
        width,
        height,
        depth);

    const StressSystem system = BuildStructuredImplicitFemLikeSystem(width, height, depth);
    const double currentMs = SolveAndTime(
        "CurrentLinearSolver",
        PhysiK::LinearSolverBackend::Current,
        system.nodes,
        system.solverData);

    if (currentMs < 0.0)
    {
        return 1;
    }

    if (PhysiK::IsLinearSolverBackendAvailable(PhysiK::LinearSolverBackend::MKL))
    {
        const double mklMs = SolveAndTime(
            "MKLLinearSolver",
            PhysiK::LinearSolverBackend::MKL,
            system.nodes,
            system.solverData);
        if (mklMs < 0.0)
        {
            return 1;
        }

        if (mklMs > 0.0)
        {
            std::printf("Speed ratio Current/MKL: %.3f\n", currentMs / mklMs);
        }
    }
    else
    {
        std::printf("MKLLinearSolver unavailable in this build; skipping MKL run.\n");
    }

    return 0;
}
