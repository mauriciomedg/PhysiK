#include "PhysiK/Geometry/TetMeshGenerator.h"

#include <cmath>
#include <cstddef>
#include <functional>
#include <unordered_map>
#include <utility>

#include "PhysiK/Math/Vec3.h"

namespace PhysiK
{
    namespace
    {
        struct CellKey
        {
            int x = 0;
            int y = 0;
            int z = 0;

            bool operator==(const CellKey& other) const
            {
                return x == other.x && y == other.y && z == other.z;
            }
        };

        struct CellKeyHash
        {
            std::size_t operator()(const CellKey& key) const
            {
                const std::size_t h0 = std::hash<int>{}(key.x);
                const std::size_t h1 = std::hash<int>{}(key.y);
                const std::size_t h2 = std::hash<int>{}(key.z);
                return h0 ^ (h1 << 1u) ^ (h2 << 2u);
            }
        };

        struct FaceKey
        {
            int node0 = -1;
            int node1 = -1;
            int node2 = -1;

            bool operator==(const FaceKey& other) const
            {
                return node0 == other.node0 &&
                    node1 == other.node1 &&
                    node2 == other.node2;
            }
        };

        struct FaceKeyHash
        {
            std::size_t operator()(const FaceKey& key) const
            {
                const std::size_t h0 = std::hash<int>{}(key.node0);
                const std::size_t h1 = std::hash<int>{}(key.node1);
                const std::size_t h2 = std::hash<int>{}(key.node2);
                return h0 ^ (h1 << 1u) ^ (h2 << 2u);
            }
        };

        CellKey MakeCellKey(const Vec3& position, float cellSize)
        {
            return CellKey{
                static_cast<int>(std::floor(position.x / cellSize)),
                static_cast<int>(std::floor(position.y / cellSize)),
                static_cast<int>(std::floor(position.z / cellSize))};
        }

        FaceKey MakeFaceKey(int node0, int node1, int node2)
        {
            FaceKey key{node0, node1, node2};
            if (key.node1 < key.node0)
            {
                std::swap(key.node0, key.node1);
            }
            if (key.node2 < key.node1)
            {
                std::swap(key.node1, key.node2);
            }
            if (key.node1 < key.node0)
            {
                std::swap(key.node0, key.node1);
            }

            return key;
        }

        bool HasDuplicateNode(int a, int b, int c, int d)
        {
            return a == b || a == c || a == d ||
                b == c || b == d ||
                c == d;
        }

        bool IsNearZeroVolume(
            const std::vector<Vec3>& positions,
            int a,
            int b,
            int c,
            int d)
        {
            const Vec3 ab = positions[static_cast<std::size_t>(b)] -
                positions[static_cast<std::size_t>(a)];
            const Vec3 ac = positions[static_cast<std::size_t>(c)] -
                positions[static_cast<std::size_t>(a)];
            const Vec3 ad = positions[static_cast<std::size_t>(d)] -
                positions[static_cast<std::size_t>(a)];
            const float sixVolume = Dot(Cross(ab, ac), ad);
            return sixVolume * sixVolume <= 1.0e-24f;
        }

        TetMeshTopologyDiagnostics BuildTopologyDiagnostics(
            const std::vector<int>& tetLocalNodeIndices)
        {
            std::unordered_map<FaceKey, int, FaceKeyHash> faceCounts;
            const int tetCount =
                static_cast<int>(tetLocalNodeIndices.size() / 4u);
            for (int tetIndex = 0; tetIndex < tetCount; ++tetIndex)
            {
                const int a = tetLocalNodeIndices[tetIndex * 4 + 0];
                const int b = tetLocalNodeIndices[tetIndex * 4 + 1];
                const int c = tetLocalNodeIndices[tetIndex * 4 + 2];
                const int d = tetLocalNodeIndices[tetIndex * 4 + 3];
                ++faceCounts[MakeFaceKey(a, b, c)];
                ++faceCounts[MakeFaceKey(a, d, b)];
                ++faceCounts[MakeFaceKey(a, c, d)];
                ++faceCounts[MakeFaceKey(b, d, c)];
            }

            TetMeshTopologyDiagnostics diagnostics;
            for (const auto& item : faceCounts)
            {
                if (item.second == 1)
                {
                    ++diagnostics.boundaryFaceCount;
                }
                else if (item.second == 2)
                {
                    ++diagnostics.internalFaceCount;
                }
                else if (item.second > 2)
                {
                    ++diagnostics.nonManifoldFaceCount;
                }
            }

            return diagnostics;
        }
    }

    GeneratedTetMesh TetMeshGenerator::Generate(
        const Vec3* positions,
        int nodeCount,
        const int* tetLocalNodeIndices,
        int tetCount,
        const TetMeshBuildOptions& options)
    {
        GeneratedTetMesh result;
        result.rawNodeCount = nodeCount > 0 ? nodeCount : 0;
        result.rawTetCount = tetCount > 0 ? tetCount : 0;
        result.weldTolerance = options.weldTolerance;

        if (positions == nullptr || nodeCount <= 0)
        {
            return result;
        }

        std::vector<int> oldNodeToNewNode(
            static_cast<std::size_t>(nodeCount),
            -1);

        const float tolerance =
            options.weldTolerance > 0.0f ? options.weldTolerance : 0.0f;
        const float toleranceSquared = tolerance * tolerance;
        const bool weldNodes =
            options.weldCoincidentNodes && tolerance > 0.0f;
        std::unordered_map<CellKey, std::vector<int>, CellKeyHash> grid;

        for (int oldNode = 0; oldNode < nodeCount; ++oldNode)
        {
            const Vec3& position = positions[oldNode];
            int existingNode = -1;

            if (weldNodes)
            {
                const CellKey cell = MakeCellKey(position, tolerance);
                for (int dz = -1; dz <= 1 && existingNode < 0; ++dz)
                {
                    for (int dy = -1; dy <= 1 && existingNode < 0; ++dy)
                    {
                        for (int dx = -1; dx <= 1 && existingNode < 0; ++dx)
                        {
                            const CellKey neighbor{
                                cell.x + dx,
                                cell.y + dy,
                                cell.z + dz};
                            const auto found = grid.find(neighbor);
                            if (found == grid.end())
                            {
                                continue;
                            }

                            for (int newNode : found->second)
                            {
                                const Vec3 delta =
                                    result.positions[static_cast<std::size_t>(newNode)] -
                                    position;
                                if (delta.LengthSquared() <= toleranceSquared)
                                {
                                    existingNode = newNode;
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            if (existingNode >= 0)
            {
                oldNodeToNewNode[static_cast<std::size_t>(oldNode)] =
                    existingNode;
                continue;
            }

            const int newNode = static_cast<int>(result.positions.size());
            result.positions.push_back(position);
            oldNodeToNewNode[static_cast<std::size_t>(oldNode)] =
                newNode;

            if (weldNodes)
            {
                grid[MakeCellKey(position, tolerance)].push_back(newNode);
            }
        }

        result.weldedNodeCount = static_cast<int>(result.positions.size());
        result.weldedAwayNodeCount = result.rawNodeCount - result.weldedNodeCount;

        if (tetLocalNodeIndices != nullptr && tetCount > 0)
        {
            result.tetLocalNodeIndices.reserve(static_cast<std::size_t>(tetCount) * 4u);
            for (int tetIndex = 0; tetIndex < tetCount; ++tetIndex)
            {
                const int raw0 = tetLocalNodeIndices[tetIndex * 4 + 0];
                const int raw1 = tetLocalNodeIndices[tetIndex * 4 + 1];
                const int raw2 = tetLocalNodeIndices[tetIndex * 4 + 2];
                const int raw3 = tetLocalNodeIndices[tetIndex * 4 + 3];
                if (raw0 < 0 || raw0 >= nodeCount ||
                    raw1 < 0 || raw1 >= nodeCount ||
                    raw2 < 0 || raw2 >= nodeCount ||
                    raw3 < 0 || raw3 >= nodeCount)
                {
                    ++result.removedDegenerateTetCount;
                    continue;
                }

                const int a = oldNodeToNewNode[static_cast<std::size_t>(raw0)];
                const int b = oldNodeToNewNode[static_cast<std::size_t>(raw1)];
                const int c = oldNodeToNewNode[static_cast<std::size_t>(raw2)];
                const int d = oldNodeToNewNode[static_cast<std::size_t>(raw3)];
                if ((options.removeDegenerateTets && HasDuplicateNode(a, b, c, d)) ||
                    (options.removeDegenerateTets &&
                        IsNearZeroVolume(result.positions, a, b, c, d)))
                {
                    ++result.removedDegenerateTetCount;
                    continue;
                }

                result.tetLocalNodeIndices.push_back(a);
                result.tetLocalNodeIndices.push_back(b);
                result.tetLocalNodeIndices.push_back(c);
                result.tetLocalNodeIndices.push_back(d);
            }
        }

        result.topologyDiagnostics =
            BuildTopologyDiagnostics(result.tetLocalNodeIndices);
        return result;
    }
}
