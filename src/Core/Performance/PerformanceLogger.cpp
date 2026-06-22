#include "PhysiK/Core/Performance/PerformanceLogger.h"

#include <filesystem>
#include <system_error>

namespace PhysiK
{
    namespace
    {
        constexpr std::size_t CgProfileFlushThreshold = 64u;
    }

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
            << record.assembleMassesMs << ','
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
            << record.assembleMassesTotalMs << ','
            << record.computeDeformationGradientMs << ','
            << record.extractRotationPolarMs << ','
            << record.averageExtractRotationPolarMs << ','
            << record.rotateElementStiffnessMs << ','
            << record.computeElasticForcesMs << ','
            << record.tetMatrixWriteMs << ','
            << record.polarCallCount << ','
            << record.averagePolarIterations << ','
            << record.maxPolarIterationsObserved << ','
            << record.polarEarlyExitCount << ','
            << record.conjugateGradientSolveMs << ','
            << record.sparseMultiplyMs << ','
            << record.preconditionerBuildMs << ','
            << record.cgTotalMs << ','
            << record.cgMultiplyMs << ','
            << record.cgApplyPreconditionerMs << ','
            << record.cgDotVectorOpsMs << ','
            << record.cgIterations << ','
            << record.cgResidual << ','
            << record.cgConverged << ','
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
                << "assembleComponentsMs,assembleMassesMs,addGravityForcesMs,"
                << "assembleConnectionsMs,"
                << "assembleComponentCount,assembleTetCount,assembleComponentTotalMs,"
                << "assembleLinearFemMs,assembleCorotationalFemMs,assembleTetForceMs,"
                << "assembleTetStiffnessMs,assembleMatrixAddBlockMs,assembleRhsWriteMs,"
                << "assembleMassesTotalMs,"
                << "computeDeformationGradientMs,extractRotationPolarMs,"
                << "averageExtractRotationPolarMs,rotateElementStiffnessMs,"
                << "computeElasticForcesMs,tetMatrixWriteMs,polarCallCount,"
                << "averagePolarIterations,maxPolarIterationsObserved,polarEarlyExitCount,"
                << "conjugateGradientSolveMs,sparseMultiplyMs,"
                << "preconditionerBuildMs,cgTotalMs,cgMultiplyMs,"
                << "cgApplyPreconditionerMs,cgDotVectorOpsMs,"
                << "cgIterations,cgResidual,cgConverged,"
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

    CgProfileCsvLogger::~CgProfileCsvLogger()
    {
        Flush();
    }

    void CgProfileCsvLogger::Log(const CgProfileRecord& record)
    {
        if (!EnsureOpen())
        {
            return;
        }

        pendingRecords.push_back(record);
        FlushIfNeeded();
    }

    bool CgProfileCsvLogger::EnsureOpen()
    {
        if (file.is_open())
        {
            return true;
        }

        constexpr const char* ProfilePath = "PhysiK_CG_Profile.csv";
        std::error_code error;
        const bool writeHeader =
            !std::filesystem::exists(ProfilePath, error) ||
            std::filesystem::file_size(ProfilePath, error) == 0u;

        file.open(ProfilePath, std::ios::out | std::ios::app);
        if (!file.is_open())
        {
            return false;
        }

        if (writeHeader)
        {
            file
                << "solve,blocks,nonZeroBlocks,maxIterations,iterations,"
                << "converged,residualNorm,tolerance,preconditionerBuildMs,"
                << "cgTotalMs,cgMultiplyMs,cgApplyPreconditionerMs,"
                << "cgDotVectorOpsMs\n";
        }

        headerWritten = true;
        return true;
    }

    void CgProfileCsvLogger::FlushIfNeeded()
    {
        if (pendingRecords.size() >= CgProfileFlushThreshold)
        {
            Flush();
        }
    }

    void CgProfileCsvLogger::Flush()
    {
        if (!file.is_open())
        {
            pendingRecords.clear();
            return;
        }

        for (const CgProfileRecord& record : pendingRecords)
        {
            ++solveCount;
            file
                << solveCount << ','
                << record.blockCount << ','
                << record.nonZeroBlockCount << ','
                << record.maxIterations << ','
                << record.iterations << ','
                << (record.converged ? 1 : 0) << ','
                << record.residualNorm << ','
                << record.tolerance << ','
                << record.preconditionerBuildMs << ','
                << record.cgTotalMs << ','
                << record.cgMultiplyMs << ','
                << record.cgApplyPreconditionerMs << ','
                << record.cgDotVectorOpsMs
                << '\n';
        }

        pendingRecords.clear();
        file.flush();
    }

    CgProfileCsvLogger& GetCgProfileCsvLogger()
    {
        static CgProfileCsvLogger logger;
        return logger;
    }
}
