#include "PhysiK/Components/TetMeshComponent.h"

#include "PhysiK/Core/Solvers/SolverData.h"
#include "PhysiK/Core/World/World.h"

namespace PhysiK
{
    void TetMeshComponent::UpdateSystem(
        World& world,
        SolverData& solverData,
        float dt)
    {
        femModel.UpdateSystem(world, *this, solverData, dt);
    }
}
