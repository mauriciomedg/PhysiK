#include "PhysiK/Core/Solvers/Linear/MKLLinearSolver.h"
#include "PhysiK/Core/Solvers/SolverData.h"
#include "PhysiK/PhysicsData/Node.h"

#if defined(PHYSIK_ENABLE_MKL)

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace
{
    void Require(bool condition, const char* message)
    {
        if (!condition)
        {
            std::fprintf(stderr, "MKLLinearSolverTests failed: %s\n", message);
            std::exit(1);
        }
    }

    bool NearlyEqual(double a, double b, double tolerance = 1.0e-9)
    {
        return std::abs(a - b) <= tolerance;
    }

    void AssertVectorNear(
        const std::vector<double>& actual,
        const std::vector<double>& expected,
        double tolerance = 1.0e-9)
    {
        Require(actual.size() == expected.size(), "double vector size mismatch");
        for (std::size_t i = 0; i < expected.size(); ++i)
        {
            Require(NearlyEqual(actual[i], expected[i], tolerance), "double vector value mismatch");
        }
    }

    PhysiK::LinearSolveSettings MakeTinySolveSettings()
    {
        PhysiK::LinearSolveSettings settings;
        settings.maxIterations = 1;
        settings.tolerance = 1.0e-10f;
        settings.useJacobiPreconditioner = false;
        return settings;
    }

    PhysiK::Mat3 DiagonalBlock(float value)
    {
        return PhysiK::Mat3::FromColumns(
            PhysiK::Vec3{value, 0.0f, 0.0f},
            PhysiK::Vec3{0.0f, value, 0.0f},
            PhysiK::Vec3{0.0f, 0.0f, value});
    }

    void AssertVectorNear(
        const std::vector<float>& actual,
        const std::vector<float>& expected,
        float tolerance = 1.0e-5f)
    {
        Require(actual.size() == expected.size(), "float vector size mismatch");
        for (std::size_t i = 0; i < expected.size(); ++i)
        {
            Require(std::abs(actual[i] - expected[i]) <= tolerance, "float vector value mismatch");
        }
    }
}

void MKLLinearSolverSolvesDiagonal2x2System()
{
    PhysiK::CSRMatrix matrix;
    matrix.rowCount = 2;
    matrix.colCount = 2;
    matrix.rowOffsets = {0, 1, 2};
    matrix.columnIndices = {0, 1};
    matrix.values = {2.0, 4.0};

    const std::vector<double> rhs = {6.0, 20.0};
    const std::vector<double> expected = {3.0, 5.0};
    std::vector<double> solution;

    PhysiK::MKLLinearSolver solver;
    const PhysiK::LinearSolveResult result =
        solver.SolveSPD(matrix, rhs, solution, MakeTinySolveSettings());

    Require(result.converged, "2x2 diagonal MKL solve did not converge");
    AssertVectorNear(solution, expected);
}

void MKLLinearSolverSolvesSmallSPD3x3System()
{
    PhysiK::CSRMatrix matrix;
    matrix.rowCount = 3;
    matrix.colCount = 3;
    matrix.rowOffsets = {0, 2, 5, 7};
    matrix.columnIndices = {
        0, 1,
        0, 1, 2,
        1, 2};
    matrix.values = {
        4.0, 1.0,
        1.0, 3.0, 1.0,
        1.0, 2.0};

    const std::vector<double> expected = {1.0, 2.0, 3.0};
    std::vector<double> rhs;
    matrix.Multiply(expected, rhs);

    std::vector<double> solution;
    PhysiK::MKLLinearSolver solver;
    const PhysiK::LinearSolveResult result =
        solver.SolveSPD(matrix, rhs, solution, MakeTinySolveSettings());

    Require(result.converged, "3x3 SPD MKL solve did not converge");
    AssertVectorNear(solution, expected);
}

void MKLBackendIsAvailableWhenCompiled()
{
    Require(
        PhysiK::IsLinearSolverBackendAvailable(PhysiK::LinearSolverBackend::MKL),
        "MKL backend is not available in MKL-enabled test build");
}

void MKLBackendSolvesSparseBlockSystemWhenSelected()
{
    PhysiK::SparseBlockMatrix matrix;
    matrix.BuildPattern(1, {{0, 0}});
    matrix.AddBlock(
        0,
        0,
        PhysiK::Mat3::FromColumns(
            PhysiK::Vec3{2.0f, 0.0f, 0.0f},
            PhysiK::Vec3{0.0f, 3.0f, 0.0f},
            PhysiK::Vec3{0.0f, 0.0f, 4.0f}));

    const std::vector<float> rhs = {2.0f, 6.0f, 12.0f};
    std::vector<float> solution;

    PhysiK::LinearSolveSettings settings;
    settings.maxIterations = 1;
    settings.tolerance = 1.0e-6f;
    settings.useJacobiPreconditioner = false;

    const PhysiK::LinearSolveResult result =
        PhysiK::GetLinearSolver(PhysiK::LinearSolverBackend::MKL).Solve(
            matrix,
            rhs,
            solution,
            settings);

    Require(result.converged, "sparse block MKL solve did not converge");
    Require(solution.size() == rhs.size(), "sparse block MKL solution size mismatch");
    Require(NearlyEqual(solution[0], 1.0, 0.0001), "sparse block MKL x value mismatch");
    Require(NearlyEqual(solution[1], 2.0, 0.0001), "sparse block MKL y value mismatch");
    Require(NearlyEqual(solution[2], 3.0, 0.0001), "sparse block MKL z value mismatch");
}

void MKLBackendSolvesPreparedSolverDataImplicitSystem()
{
    std::vector<PhysiK::Node> nodes(2);
    nodes[0].fixed = false;
    nodes[1].fixed = false;

    PhysiK::SolverData solverData;
    solverData.AddNodeMass(0, 2.0f);
    solverData.AddNodeMass(1, 3.0f);
    solverData.AddNodeForce(0, PhysiK::Vec3{8.0f, 16.0f, 24.0f});
    solverData.AddNodeForce(1, PhysiK::Vec3{24.0f, 30.0f, 36.0f});
    solverData.AddStiffnessBlock(0, 0, DiagonalBlock(200.0f));
    solverData.AddStiffnessBlock(0, 1, DiagonalBlock(-50.0f));
    solverData.AddStiffnessBlock(1, 0, DiagonalBlock(-50.0f));
    solverData.AddStiffnessBlock(1, 1, DiagonalBlock(100.0f));

    constexpr float Dt = 0.1f;
    Require(
        solverData.GetLinearSolverBackend() == PhysiK::LinearSolverBackend::Current,
        "SolverData default backend is not Current");
    Require(solverData.PrecomputeImplicitSolve(nodes, Dt), "SolverData precompute failed");

    solverData.SetLinearSolverBackend(PhysiK::LinearSolverBackend::MKL);
    Require(
        solverData.GetLinearSolverBackend() == PhysiK::LinearSolverBackend::MKL,
        "SolverData did not store MKL backend selection");
    Require(
        solverData.SolveImplicitLinearSystem(),
        "SolverData implicit linear solve with MKL failed");

    const std::vector<float> expected = {
        0.279365f,
        0.501587f,
        0.723810f,
        0.634921f,
        0.812698f,
        0.990476f};
    AssertVectorNear(solverData.GetDeltaVelocity(), expected);
}

int main()
{
    MKLLinearSolverSolvesDiagonal2x2System();
    MKLLinearSolverSolvesSmallSPD3x3System();
    MKLBackendIsAvailableWhenCompiled();
    MKLBackendSolvesSparseBlockSystemWhenSelected();
    MKLBackendSolvesPreparedSolverDataImplicitSystem();
    return 0;
}

#else

int main()
{
    return 0;
}

#endif
