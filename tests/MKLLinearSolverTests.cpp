#include "PhysiK/Core/Solvers/Linear/MKLLinearSolver.h"

#if defined(PHYSIK_ENABLE_MKL)

#include <cassert>
#include <cmath>
#include <vector>

namespace
{
    bool NearlyEqual(double a, double b, double tolerance = 1.0e-9)
    {
        return std::abs(a - b) <= tolerance;
    }

    void AssertVectorNear(
        const std::vector<double>& actual,
        const std::vector<double>& expected,
        double tolerance = 1.0e-9)
    {
        assert(actual.size() == expected.size());
        for (std::size_t i = 0; i < expected.size(); ++i)
        {
            assert(NearlyEqual(actual[i], expected[i], tolerance));
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

    assert(result.converged);
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

    assert(result.converged);
    AssertVectorNear(solution, expected);
}

int main()
{
    MKLLinearSolverSolvesDiagonal2x2System();
    MKLLinearSolverSolvesSmallSPD3x3System();
    return 0;
}

#else

int main()
{
    return 0;
}

#endif
