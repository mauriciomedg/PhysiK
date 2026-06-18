#include "PhysiK/Core/Solvers/Linear/ConjugateGradientSolver.h"

#include <algorithm>
#include <cmath>

namespace PhysiK
{
    namespace
    {
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

        void ApplyPreconditioner(
            const std::vector<Mat3>& MInv,
            const std::vector<Vec3>& input,
            std::vector<Vec3>& output)
        {
            output.assign(input.size(), Vec3{});

            if (MInv.size() != input.size())
            {
                output = input;
                return;
            }

            for (std::size_t i = 0; i < input.size(); ++i)
            {
                output[i] = MInv[i] * input[i];
            }
        }

        float ResidualNorm(const std::vector<Vec3>& residual)
        {
            const float squared = Dot(residual, residual);
            if (!IsFinite(squared) || squared < 0.0f)
            {
                return 0.0f;
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
        ApplyPreconditioner(MInv, r, d);

        if (!IsFinite(d))
        {
            return result;
        }

        float deltaNew = Dot(r, d);
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
            result.residualNorm = ResidualNorm(r);
            return result;
        }

        for (int iteration = 0; iteration < maxIterations; ++iteration)
        {
            if (deltaNew <= target)
            {
                result.converged = true;
                break;
            }

            A.Multiply(d, qOrS);
            if (qOrS.size() != b.size() || !IsFinite(qOrS))
            {
                break;
            }

            const float dDotQ = Dot(d, qOrS);
            if (!IsFinite(dDotQ) || dDotQ <= 0.0f)
            {
                break;
            }

            const float alpha = deltaNew / dDotQ;
            if (!IsFinite(alpha))
            {
                break;
            }

            AddScaled(x, d, alpha);
            AddScaled(r, qOrS, -alpha);

            ApplyPreconditioner(MInv, r, qOrS);

            const float deltaOld = deltaNew;
            deltaNew = Dot(r, qOrS);
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

            for (std::size_t i = 0; i < d.size(); ++i)
            {
                d[i] = qOrS[i] + d[i] * beta;
            }

            result.iterations = iteration + 1;
        }

        if (deltaNew <= target)
        {
            result.converged = true;
        }

        result.residualNorm = ResidualNorm(r);
        return result;
    }
}
