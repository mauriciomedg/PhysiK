#include "PhysiK/Core/Performance/PerformanceLogger.h"

#include <filesystem>
#include <system_error>

namespace PhysiK
{
    PerformanceTimer::PerformanceTimer()
        : start(Clock::now())
    {
    }

    double PerformanceTimer::ElapsedMilliseconds() const
    {
        return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
    }

    void PerformanceLogger::Enable(bool value)
    {
        enabled = value;
        if (!enabled)
        {
            Close();
        }
    }

    bool PerformanceLogger::IsEnabled() const
    {
        return enabled;
    }

    void PerformanceLogger::SetPath(const char* value)
    {
        path = value != nullptr && value[0] != '\0' ?
            std::string(value) :
            std::string("logs/physik_performance.csv");
        Close();
    }

    const std::string& PerformanceLogger::GetPath() const
    {
        return path;
    }

    void PerformanceLogger::Log(const PerformanceLogRecord& record)
    {
        if (!enabled || !EnsureOpen())
        {
            return;
        }

        file
            << record.frameIndex << ','
            << record.substepIndex << ','
            << record.dt << ','
            << record.totalStepMs << ','
            << record.generateCollisionConnectionsMs << ','
            << record.buildSolverDataMs << ','
            << record.assembleSystemMs << ','
            << record.assembleComponentsMs << ','
            << record.addDefaultMassesMs << ','
            << record.addGravityForcesMs << ','
            << record.assembleConnectionsMs << ','
            << record.assembleComponentCount << ','
            << record.assembleTetCount << ','
            << record.assembleComponentTotalMs << ','
            << record.assembleLinearFemMs << ','
            << record.assembleCorotationalFemMs << ','
            << record.assembleTetForceMs << ','
            << record.assembleTetStiffnessMs << ','
            << record.assembleMatrixAddBlockMs << ','
            << record.assembleRhsWriteMs << ','
            << record.addDefaultMassesTotalMs << ','
            << record.addDefaultMassLoopMs << ','
            << record.addDefaultMassMatrixWriteMs << ','
            << record.addDefaultMassDynamicNodeCount << ','
            << record.conjugateGradientSolveMs << ','
            << record.sparseMultiplyMs << ','
            << record.cgIterations << ','
            << record.cgResidual << ','
            << record.dynamicBlockCount << ','
            << record.tetCount << ','
            << record.activeTetCount << ','
            << record.transientConnectionCount
            << '\n';
        file.flush();
    }

    bool PerformanceLogger::EnsureOpen()
    {
        if (file.is_open())
        {
            return true;
        }

        std::error_code error;
        const std::filesystem::path logPath(path);
        const std::filesystem::path parentPath = logPath.parent_path();
        if (!parentPath.empty())
        {
            std::filesystem::create_directories(parentPath, error);
        }

        const bool writeHeader =
            !std::filesystem::exists(logPath, error) ||
            std::filesystem::file_size(logPath, error) == 0u;
        file.open(path, std::ios::out | std::ios::app);
        if (!file.is_open())
        {
            return false;
        }

        if (writeHeader)
        {
            file
                << "frameIndex,substepIndex,dt,totalStepMs,"
                << "generateCollisionConnectionsMs,buildSolverDataMs,assembleSystemMs,"
                << "assembleComponentsMs,addDefaultMassesMs,addGravityForcesMs,"
                << "assembleConnectionsMs,"
                << "assembleComponentCount,assembleTetCount,assembleComponentTotalMs,"
                << "assembleLinearFemMs,assembleCorotationalFemMs,assembleTetForceMs,"
                << "assembleTetStiffnessMs,assembleMatrixAddBlockMs,assembleRhsWriteMs,"
                << "addDefaultMassesTotalMs,addDefaultMassLoopMs,"
                << "addDefaultMassMatrixWriteMs,addDefaultMassDynamicNodeCount,"
                << "conjugateGradientSolveMs,sparseMultiplyMs,cgIterations,cgResidual,"
                << "dynamicBlockCount,tetCount,activeTetCount,transientConnectionCount\n";
        }
        headerWritten = true;

        return true;
    }

    void PerformanceLogger::Close()
    {
        if (file.is_open())
        {
            file.close();
        }
        headerWritten = false;
    }
}
