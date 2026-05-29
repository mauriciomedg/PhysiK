#pragma once

#include <vector>

#include "PhysiK/API/PhysiKAPI.h"
#include "PhysiK/Math/Vec3.h"

namespace PhysiK
{
    struct TetMeshBuildOptions
    {
        bool weldCoincidentNodes = true;
        float weldTolerance = 1.0e-5f;
        bool removeDegenerateTets = true;
    };

    struct TetMeshTopologyDiagnostics
    {
        int boundaryFaceCount = 0;
        int internalFaceCount = 0;
        int nonManifoldFaceCount = 0;
    };

    struct GeneratedTetMesh
    {
        std::vector<Vec3> positions;
        std::vector<int> tetLocalNodeIndices;
        TetMeshTopologyDiagnostics topologyDiagnostics;
        int rawNodeCount = 0;
        int rawTetCount = 0;
        int weldedNodeCount = 0;
        int weldedAwayNodeCount = 0;
        int removedDegenerateTetCount = 0;
        float weldTolerance = 0.0f;
    };

    class PHYSIK_API TetMeshGenerator
    {
    public:
        static GeneratedTetMesh Generate(
            const Vec3* positions,
            int nodeCount,
            const int* tetLocalNodeIndices,
            int tetCount,
            const TetMeshBuildOptions& options = TetMeshBuildOptions{});
    };
}
