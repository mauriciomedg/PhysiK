#include "PhysiK/API/PhysiKAPI.h"

#include <cassert>
#include <cmath>

namespace
{
    struct Point
    {
        float x = 0.0f;
        float y = 0.0f;
        float z = 0.0f;
    };

    Point GetNodePosition(PhysiK::WorldHandle world, int nodeIndex)
    {
        Point point;
        PHYSIK_GetNodePosition(world, nodeIndex, &point.x, &point.y, &point.z);
        return point;
    }

    Point BarycentricPoint(
        const Point& p0,
        const Point& p1,
        const Point& p2,
        const Point& p3,
        float w0,
        float w1,
        float w2,
        float w3)
    {
        return Point{
            p0.x * w0 + p1.x * w1 + p2.x * w2 + p3.x * w3,
            p0.y * w0 + p1.y * w1 + p2.y * w2 + p3.y * w3,
            p0.z * w0 + p1.z * w1 + p2.z * w2 + p3.z * w3};
    }

    float DistanceSquared(const Point& a, const Point& b)
    {
        const float dx = a.x - b.x;
        const float dy = a.y - b.y;
        const float dz = a.z - b.z;
        return dx * dx + dy * dy + dz * dz;
    }
}

int main()
{
    PhysiK::WorldHandle world = PHYSIK_CreateWorld();
    assert(world != nullptr);

    const int node0 = PHYSIK_AddNode(world, 0.0f, 0.0f, 0.0f, 1.0f);
    const int node1 = PHYSIK_AddNode(world, 1.0f, 0.0f, 0.0f, 1.0f);
    const int node2 = PHYSIK_AddNode(world, 0.0f, 1.0f, 0.0f, 1.0f);
    const int node3 = PHYSIK_AddNode(world, 0.0f, 0.0f, 1.0f, 1.0f);
    PHYSIK_AddTet(world, node0, node1, node2, node3);

    const int nodes[] = {node0, node1, node2, node3};
    const int tets[] = {0};
    assert(PHYSIK_CreateTetMeshComponent(world, nodes, 4, tets, 1) != nullptr);

    const float w0 = 0.25f;
    const float w1 = 0.25f;
    const float w2 = 0.25f;
    const float w3 = 0.25f;
    const Point target{0.25f, 0.25f, 1.25f};

    const Point before = BarycentricPoint(
        GetNodePosition(world, node0),
        GetNodePosition(world, node1),
        GetNodePosition(world, node2),
        GetNodePosition(world, node3),
        w0,
        w1,
        w2,
        w3);

    PHYSIK_AddPointConnection(
        world,
        node0,
        node1,
        node2,
        node3,
        w0,
        w1,
        w2,
        w3,
        target.x,
        target.y,
        target.z,
        100.0f,
        0.0f);

    PHYSIK_Step(world, 0.1f);

    const Point after = BarycentricPoint(
        GetNodePosition(world, node0),
        GetNodePosition(world, node1),
        GetNodePosition(world, node2),
        GetNodePosition(world, node3),
        w0,
        w1,
        w2,
        w3);

    assert(DistanceSquared(after, target) < DistanceSquared(before, target));
    assert(after.z > before.z);

    PHYSIK_DestroyWorld(world);
    return 0;
}
