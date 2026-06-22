#include "PhysiK/Core/Solvers/Linear/ConjugateGradientSolver.h"

#include <algorithm>
#include <cmath>
#include <limits>

#if defined(PHYSIK_ENABLE_PERF_LOGGING) || defined(PHYSIK_ENABLE_SOLVER_PROFILING)
#include <chrono>
#define PHYSIK_COLLECT_CG_TIMING 1
#endif

namespace PhysiK
{
    namespace
    {
#if defined(PHYSIK_COLLECT_CG_TIMING)
        using Clock = std::chrono::steady_clock;

        double ElapsedMilliseconds(Clock::time_point start)
        {
            return std::chrono::duration<double, std::milli>(
                Clock::now() - start).count();
        }
#endif

        bool IsFinite(float value)
        {
            return std::isfinite(value);
        }

        bool IsFinite(const Vec3& value)
        {
            return IsFinite(value.x) && IsFinite(value.y) && IsFinite(value.z);
        }

        bool IsFinite(const std::vector<Vec3>& values)
        {
            for (const Vec3& value : values)
            {
                if (!IsFinite(value))
                {
                    return false;
                }
            }

            return true;
        }

        bool IsFinite(const Mat3& matrix)
        {
            return IsFinite(matrix.columns[0]) &&
                IsFinite(matrix.columns[1]) &&
                IsFinite(matrix.columns[2]);
        }

        float Dot(const std::vector<Vec3>& a, const std::vector<Vec3>& b)
        {
            float result = 0.0f;
            const std::size_t count = std::min(a.size(), b.size());

            for (std::size_t i = 0; i < count; ++i)
            {
                result += PhysiK::Dot(a[i], b[i]);
            }

            return result;
        }

        void AddScaled(
            std::vector<Vec3>& dst,
            const std::vector<Vec3>& src,
            float scale)
        {
            const std::size_t count = std::min(dst.size(), src.size());

            for (std::size_t i = 0; i < count; ++i)
            {
                dst[i] += src[i] * scale;
            }
        }

        bool ApplyPreconditioner(
            const std::vector<Mat3>& MInv,
            const std::vector<Vec3>& input,
            std::vector<Vec3>& output)
        {
            if (MInv.size() != input.size())
            {
                output.clear();
                return false;
            }

            if (output.size() != input.size())
            {
                output.resize(input.size());
            }

            for (std::size_t i = 0; i < input.size(); ++i)
            {
                if (!IsFinite(MInv[i]))
                {
                    output.clear();
                    return false;
                }

                output[i] = MInv[i] * input[i];
                if (!IsFinite(output[i]))
                {
                    output.clear();
                    return false;
                }
            }

            return true;
        }

        float ResidualNorm(const std::vector<Vec3>& residual)
        {
            const float squared = Dot(residual, residual);
            if (!IsFinite(squared) || squared < 0.0f)
            {
                return std::numeric_limits<float>::infinity();
            }

            return std::sqrt(squared);
        }
    }

    ConjugateGradientResult SolvePreconditionedConjugateGradient(
        std::vector<Vec3>& x,
        const SparseBlockMatrix& A,
        const std::vector<Vec3>& b,
        int maxIterations,
        float tolerance,
        const std::vector<Mat3>& MInv,
        std::vector<Vec3>& r,
        std::vector<Vec3>& d,
        std::vector<Vec3>& qOrS)
    {
        ConjugateGradientResult result;

        if (b.empty())
        {
            x.clear();
            r.clear();
            d.clear();
            qOrS.clear();
            result.converged = true;
            return result;
        }

        if (A.blockCount != static_cast<int>(b.size()) ||
            A.rowStart.size() != static_cast<std::size_t>(A.blockCount + 1) ||
            maxIterations <= 0 ||
            !IsFinite(tolerance) ||
            tolerance <= 0.0f ||
            !IsFinite(b))
        {
            x.clear();
            r.clear();
            d.clear();
            qOrS.clear();
            return result;
        }

        x.assign(b.size(), Vec3{});
        r = b;
#if defined(PHYSIK_COLLECT_CG_TIMING)
        Clock::time_point timerStart = Clock::now();
#endif
        if (!ApplyPreconditioner(MInv, r, d))
        {
#if defined(PHYSIK_COLLECT_CG_TIMING)
            result.cgApplyPreconditionerMs += ElapsedMilliseconds(timerStart);
#endif
            return result;
        }
#if defined(PHYSIK_COLLECT_CG_TIMING)
        result.cgApplyPreconditionerMs += ElapsedMilliseconds(timerStart);

        timerStart = Clock::now();
#endif
        float deltaNew = Dot(r, d);
#if defined(PHYSIK_COLLECT_CG_TIMING)
        result.cgDotVectorOpsMs += ElapsedMilliseconds(timerStart);
#endif
        if (!IsFinite(deltaNew) || deltaNew < 0.0f)
        {
            return result;
        }

        const float initialDelta = deltaNew;
        const float target = tolerance * tolerance * initialDelta;

        if (!IsFinite(target))
        {
            return result;
        }

        if (deltaNew <= target)
        {
            result.converged = true;
#if defined(PHYSIK_COLLECT_CG_TIMING)
            timerStart = Clock::now();
#endif
            result.residualNorm = ResidualNorm(r);
#if defined(PHYSIK_COLLECT_CG_TIMING)
            result.cgDotVectorOpsMs += ElapsedMilliseconds(timerStart);
#endif
            return result;
        }

        for (int iteration = 0; iteration < maxIterations; ++iteration)
        {
            if (deltaNew <= target)
            {
                result.converged = true;
                break;
            }

#if defined(PHYSIK_COLLECT_CG_TIMING)
            timerStart = Clock::now();
#endif
            A.Multiply(d, qOrS);
#if defined(PHYSIK_COLLECT_CG_TIMING)
            result.cgMultiplyMs += ElapsedMilliseconds(timerStart);
#endif
            if (qOrS.size() != b.size() || !IsFinite(qOrS))
            {
                break;
            }

#if defined(PHYSIK_COLLECT_CG_TIMING)
            timerStart = Clock::now();
#endif
            const float dDotQ = Dot(d, qOrS);
#if defined(PHYSIK_COLLECT_CG_TIMING)
            result.cgDotVectorOpsMs += ElapsedMilliseconds(timerStart);
#endif
            if (!IsFinite(dDotQ) || dDotQ <= 0.0f)
            {
                break;
            }

            const float alpha = deltaNew / dDotQ;
            if (!IsFinite(alpha))
            {
                break;
            }

#if defined(PHYSIK_COLLECT_CG_TIMING)
            timerStart = Clock::now();
#endif
            AddScaled(x, d, alpha);
            AddScaled(r, qOrS, -alpha);
#if defined(PHYSIK_COLLECT_CG_TIMING)
            result.cgDotVectorOpsMs += ElapsedMilliseconds(timerStart);

            timerStart = Clock::now();
#endif
            if (!ApplyPreconditioner(MInv, r, qOrS))
            {
#if defined(PHYSIK_COLLECT_CG_TIMING)
                result.cgApplyPreconditionerMs += ElapsedMilliseconds(timerStart);
#endif
                break;
            }
#if defined(PHYSIK_COLLECT_CG_TIMING)
            result.cgApplyPreconditionerMs += ElapsedMilliseconds(timerStart);
#endif

            const float deltaOld = deltaNew;
#if defined(PHYSIK_COLLECT_CG_TIMING)
            timerStart = Clock::now();
#endif
            deltaNew = Dot(r, qOrS);
#if defined(PHYSIK_COLLECT_CG_TIMING)
            result.cgDotVectorOpsMs += ElapsedMilliseconds(timerStart);
#endif
            if (!IsFinite(deltaNew) ||
                !IsFinite(deltaOld) ||
                deltaOld <= 0.0f)
            {
                break;
            }

            const float beta = deltaNew / deltaOld;
            if (!IsFinite(beta))
            {
                break;
            }

#if defined(PHYSIK_COLLECT_CG_TIMING)
            timerStart = Clock::now();
#endif
            for (std::size_t i = 0; i < d.size(); ++i)
            {
                d[i] = qOrS[i] + d[i] * beta;
            }
#if defined(PHYSIK_COLLECT_CG_TIMING)
            result.cgDotVectorOpsMs += ElapsedMilliseconds(timerStart);
#endif

            result.iterations = iteration + 1;
        }

        if (deltaNew <= target)
        {
            result.converged = true;
        }

#if defined(PHYSIK_COLLECT_CG_TIMING)
        timerStart = Clock::now();
#endif
        result.residualNorm = ResidualNorm(r);
#if defined(PHYSIK_COLLECT_CG_TIMING)
        result.cgDotVectorOpsMs += ElapsedMilliseconds(timerStart);
#endif
        return result;
    }
}
