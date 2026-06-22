#pragma once

#include <chrono>
#include <cstdint>
#include <fstream>
#include <string>

#include "PhysiK/API/PhysiKAPI.h"

namespace PhysiK
{
    struct PerformanceLogRecord
    {
        std::uint64_t frameIndex = 0;
        int substepIndex = 0;
        float dt = 0.0f;
        double totalStepMs = 0.0;
        double generateCollisionConnectionsMs = 0.0;
        double buildSolverDataMs = 0.0;
        double assembleSystemMs = 0.0;
        double assembleComponentsMs = 0.0;
        double assembleMassesMs = 0.0;
        double addGravityForcesMs = 0.0;
        double assembleConnectionsMs = 0.0;
        int assembleComponentCount = 0;
        int assembleTetCount = 0;
        double assembleComponentTotalMs = 0.0;
        double assembleLinearFemMs = 0.0;
        double assembleCorotationalFemMs = 0.0;
        double assembleTetForceMs = 0.0;
        double assembleTetStiffnessMs = 0.0;
        double assembleMatrixAddBlockMs = 0.0;
        double assembleRhsWriteMs = 0.0;
        double assembleMassesTotalMs = 0.0;
        double computeDeformationGradientMs = 0.0;
        double extractRotationPolarMs = 0.0;
        double averageExtractRotationPolarMs = 0.0;
        double rotateElementStiffnessMs = 0.0;
        double computeElasticForcesMs = 0.0;
        double tetMatrixWriteMs = 0.0;
        int polarCallCount = 0;
        double averagePolarIterations = 0.0;
        int maxPolarIterationsObserved = 0;
        int polarEarlyExitCount = 0;
        double conjugateGradientSolveMs = 0.0;
        double sparseMultiplyMs = 0.0;
        double preconditionerBuildMs = 0.0;
        double cgTotalMs = 0.0;
        double cgMultiplyMs = 0.0;
        double cgApplyPreconditionerMs = 0.0;
        double cgDotVectorOpsMs = 0.0;
        int cgIterations = 0;
        float cgResidual = 0.0f;
        int cgConverged = 0;
        int dynamicBlockCount = 0;
        int tetCount = 0;
        int activeTetCount = 0;
        int transientConnectionCount = 0;
        double assembleStiffnessBlocksMs = 0.0;
    };

    class PHYSIK_API PerformanceTimer
    {
    public:
        PerformanceTimer();

        double ElapsedMilliseconds() const;

    private:
        using Clock = std::chrono::steady_clock;
        Clock::time_point start;
    };

    class PHYSIK_API PerformanceLogger
    {
    public:
        void Enable(bool value);
        bool IsEnabled() const;

        void SetPath(const char* value);
        const std::string& GetPath() const;

        void Log(const PerformanceLogRecord& record);

    private:
        bool EnsureOpen();
        void Close();

        bool enabled = true;
        bool headerWritten = false;
        std::string path = "logs/physik_performance.csv";
        std::ofstream file;
    };
}
