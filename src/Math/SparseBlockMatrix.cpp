#include "PhysiK/Math/SparseBlockMatrix.h"

#include <algorithm>
#include <condition_variable>
#include <mutex>
#include <thread>

namespace PhysiK
{
    namespace
    {
        constexpr int ParallelMultiplyMinBlockRows = 8192;
        constexpr int MinimumRowsPerMultiplyWorker = 2048;
        constexpr int MaxParallelMultiplyThreads = 4;

        Mat3 Add(const Mat3& a, const Mat3& b)
        {
            return Mat3::FromColumns(
                a.columns[0] + b.columns[0],
                a.columns[1] + b.columns[1],
                a.columns[2] + b.columns[2]);
        }

        bool IsValidCoordinate(int rowBlock, int colBlock, int blockCount)
        {
            return rowBlock >= 0 && rowBlock < blockCount &&
                colBlock >= 0 && colBlock < blockCount;
        }

        int GetMultiplyThreadCount(int blockCount)
        {
            if (blockCount < ParallelMultiplyMinBlockRows)
            {
                return 1;
            }

            const unsigned int hardwareThreads = std::thread::hardware_concurrency();
            const int availableThreads =
                hardwareThreads > 0u ? static_cast<int>(hardwareThreads) : 1;
            const int usefulThreads =
                std::max(1, blockCount / MinimumRowsPerMultiplyWorker);
            return std::max(
                1,
                std::min({availableThreads, usefulThreads, MaxParallelMultiplyThreads}));
        }

        void MultiplyBlockRows(
            int rowBeginBlock,
            int rowEndBlock,
            const float* inputValues,
            float* outputValues,
            const int* rowStarts,
            const int* columnIndices,
            const Mat3* blockValues)
        {
            for (int rowBlock = rowBeginBlock; rowBlock < rowEndBlock; ++rowBlock)
            {
                float rowX = 0.0f;
                float rowY = 0.0f;
                float rowZ = 0.0f;
                const int rowBegin = rowStarts[rowBlock];
                const int rowEnd = rowStarts[rowBlock + 1];
                for (int blockIndex = rowBegin; blockIndex < rowEnd; ++blockIndex)
                {
                    const int columnBase = columnIndices[blockIndex] * 3;
                    const float x = inputValues[columnBase + 0];
                    const float y = inputValues[columnBase + 1];
                    const float z = inputValues[columnBase + 2];
                    const Mat3& block = blockValues[blockIndex];
                    const Vec3& column0 = block.columns[0];
                    const Vec3& column1 = block.columns[1];
                    const Vec3& column2 = block.columns[2];

                    rowX += column0.x * x + column1.x * y + column2.x * z;
                    rowY += column0.y * x + column1.y * y + column2.y * z;
                    rowZ += column0.z * x + column1.z * y + column2.z * z;
                }

                const int rowBase = rowBlock * 3;
                outputValues[rowBase + 0] = rowX;
                outputValues[rowBase + 1] = rowY;
                outputValues[rowBase + 2] = rowZ;
            }
        }

        class MultiplyThreadPool
        {
        public:
            explicit MultiplyThreadPool(int workerCount)
            {
                workers.reserve(static_cast<std::size_t>(std::max(0, workerCount)));
                for (int workerIndex = 0; workerIndex < workerCount; ++workerIndex)
                {
                    workers.emplace_back([this, workerIndex]() { WorkerLoop(workerIndex); });
                }
            }

            ~MultiplyThreadPool()
            {
                {
                    std::lock_guard<std::mutex> lock(mutex);
                    stopping = true;
                    ++generation;
                }
                workAvailable.notify_all();

                for (std::thread& worker : workers)
                {
                    if (worker.joinable())
                    {
                        worker.join();
                    }
                }
            }

            int WorkerCapacity() const
            {
                return static_cast<int>(workers.size());
            }

            void Run(
                int threadCount,
                int blockCount,
                const float* inputValues,
                float* outputValues,
                const int* rowStarts,
                const int* columnIndices,
                const Mat3* blockValues)
            {
                const int workerCount = std::min(threadCount - 1, WorkerCapacity());
                if (workerCount <= 0)
                {
                    MultiplyBlockRows(
                        0,
                        blockCount,
                        inputValues,
                        outputValues,
                        rowStarts,
                        columnIndices,
                        blockValues);
                    return;
                }

                const int rowsPerThread = (blockCount + threadCount - 1) / threadCount;
                {
                    std::lock_guard<std::mutex> lock(mutex);
                    taskBlockCount = blockCount;
                    taskRowsPerThread = rowsPerThread;
                    taskWorkerCount = workerCount;
                    taskInputValues = inputValues;
                    taskOutputValues = outputValues;
                    taskRowStarts = rowStarts;
                    taskColumnIndices = columnIndices;
                    taskBlockValues = blockValues;
                    remainingWorkers = workerCount;
                    ++generation;
                }
                workAvailable.notify_all();

                const int mainRowBegin = workerCount * rowsPerThread;
                MultiplyBlockRows(
                    mainRowBegin,
                    blockCount,
                    inputValues,
                    outputValues,
                    rowStarts,
                    columnIndices,
                    blockValues);

                std::unique_lock<std::mutex> lock(mutex);
                workFinished.wait(lock, [this]() { return remainingWorkers == 0; });
            }

        private:
            void WorkerLoop(int workerIndex)
            {
                int observedGeneration = 0;
                for (;;)
                {
                    int rowBegin = 0;
                    int rowEnd = 0;
                    const float* inputValues = nullptr;
                    float* outputValues = nullptr;
                    const int* rowStarts = nullptr;
                    const int* columnIndices = nullptr;
                    const Mat3* blockValues = nullptr;
                    bool hasWork = false;

                    {
                        std::unique_lock<std::mutex> lock(mutex);
                        workAvailable.wait(
                            lock,
                            [this, observedGeneration]()
                            {
                                return stopping || generation != observedGeneration;
                            });

                        if (stopping)
                        {
                            return;
                        }

                        observedGeneration = generation;
                        hasWork = workerIndex < taskWorkerCount;
                        if (hasWork)
                        {
                            rowBegin = workerIndex * taskRowsPerThread;
                            rowEnd = std::min(taskBlockCount, rowBegin + taskRowsPerThread);
                            inputValues = taskInputValues;
                            outputValues = taskOutputValues;
                            rowStarts = taskRowStarts;
                            columnIndices = taskColumnIndices;
                            blockValues = taskBlockValues;
                        }
                    }

                    if (hasWork)
                    {
                        MultiplyBlockRows(
                            rowBegin,
                            rowEnd,
                            inputValues,
                            outputValues,
                            rowStarts,
                            columnIndices,
                            blockValues);

                        std::lock_guard<std::mutex> lock(mutex);
                        --remainingWorkers;
                        if (remainingWorkers == 0)
                        {
                            workFinished.notify_one();
                        }
                    }
                }
            }

            std::vector<std::thread> workers;
            std::mutex mutex;
            std::condition_variable workAvailable;
            std::condition_variable workFinished;
            bool stopping = false;
            int generation = 0;
            int remainingWorkers = 0;
            int taskBlockCount = 0;
            int taskRowsPerThread = 0;
            int taskWorkerCount = 0;
            const float* taskInputValues = nullptr;
            float* taskOutputValues = nullptr;
            const int* taskRowStarts = nullptr;
            const int* taskColumnIndices = nullptr;
            const Mat3* taskBlockValues = nullptr;
        };

        MultiplyThreadPool& GetMultiplyThreadPool()
        {
            thread_local MultiplyThreadPool threadPool(MaxParallelMultiplyThreads - 1);
            return threadPool;
        }
    }

    void SparseBlockMatrix::Clear()
    {
        blockCount = 0;
        rowStart.clear();
        colIndex.clear();
        values.clear();
        blockLookup.clear();
    }

    void SparseBlockMatrix::ClearValues()
    {
        std::fill(values.begin(), values.end(), Mat3::Zero());
    }

    void SparseBlockMatrix::BuildPattern(
        int newBlockCount,
        const std::vector<std::pair<int, int>>& blockCoordinates)
    {
        Clear();
        if (newBlockCount <= 0)
        {
            return;
        }

        blockCount = newBlockCount;
        std::vector<std::vector<int>> rows(static_cast<std::size_t>(blockCount));
        for (const std::pair<int, int>& coordinate : blockCoordinates)
        {
            const int rowBlock = coordinate.first;
            const int colBlock = coordinate.second;
            if (!IsValidCoordinate(rowBlock, colBlock, blockCount))
            {
                continue;
            }

            rows[static_cast<std::size_t>(rowBlock)].push_back(colBlock);
        }

        rowStart.assign(static_cast<std::size_t>(blockCount + 1), 0);
        for (int row = 0; row < blockCount; ++row)
        {
            std::vector<int>& columns = rows[static_cast<std::size_t>(row)];
            std::sort(columns.begin(), columns.end());
            columns.erase(std::unique(columns.begin(), columns.end()), columns.end());
            rowStart[static_cast<std::size_t>(row + 1)] =
                rowStart[static_cast<std::size_t>(row)] + static_cast<int>(columns.size());

            for (int column : columns)
            {
                const int blockIndex = static_cast<int>(colIndex.size());
                colIndex.push_back(column);
                values.push_back(Mat3::Zero());
                blockLookup[MakeKey(row, column)] = blockIndex;
            }
        }
    }

    bool SparseBlockMatrix::AddBlock(int rowBlock, int colBlock, const Mat3& block)
    {
        const int blockIndex = FindBlockIndex(rowBlock, colBlock);
        if (blockIndex < 0)
        {
            return false;
        }

        values[static_cast<std::size_t>(blockIndex)] =
            Add(values[static_cast<std::size_t>(blockIndex)], block);
        return true;
    }

    void SparseBlockMatrix::Multiply(
        const std::vector<float>& input,
        std::vector<float>& output) const
    {
        const std::size_t dimension = static_cast<std::size_t>(std::max(0, blockCount) * 3);
        if (input.size() < dimension || rowStart.size() != static_cast<std::size_t>(blockCount + 1))
        {
            output.assign(dimension, 0.0f);
            return;
        }

        output.resize(dimension);

        const float* inputValues = input.data();
        float* outputValues = output.data();
        const int* rowStarts = rowStart.data();
        const int* columnIndices = colIndex.data();
        const Mat3* blockValues = values.data();

        const int threadCount = GetMultiplyThreadCount(blockCount);
        if (threadCount <= 1)
        {
            MultiplyBlockRows(
                0,
                blockCount,
                inputValues,
                outputValues,
                rowStarts,
                columnIndices,
                blockValues);
            return;
        }

        GetMultiplyThreadPool().Run(
            threadCount,
            blockCount,
            inputValues,
            outputValues,
            rowStarts,
            columnIndices,
            blockValues);
    }

    int SparseBlockMatrix::FindBlockIndex(int rowBlock, int colBlock) const
    {
        const auto it = blockLookup.find(MakeKey(rowBlock, colBlock));
        if (it == blockLookup.end())
        {
            return -1;
        }

        return it->second;
    }

    std::uint64_t SparseBlockMatrix::MakeKey(int rowBlock, int colBlock)
    {
        return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(rowBlock)) << 32u) |
            static_cast<std::uint32_t>(colBlock);
    }
}
