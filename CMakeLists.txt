cmake_minimum_required(VERSION 3.20)

project(PhysiK VERSION 0.1.0 LANGUAGES CXX)

option(PHYSIK_BUILD_TESTS "Build PhysiK tests" ON)
option(PHYSIK_BUILD_BENCHMARKS "Build PhysiK standalone benchmark executables" OFF)
option(PHYSIK_ENABLE_PERF_LOGGING "Enable optional CSV performance logging instrumentation" OFF)

add_library(PhysiK SHARED)
add_library(PhysiK::PhysiK ALIAS PhysiK)

target_compile_features(PhysiK PUBLIC cxx_std_17)
target_include_directories(
    PhysiK
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_sources(
    PhysiK
    PUBLIC
        include/PhysiK/API/Handles.h
        include/PhysiK/API/PhysiKAPI.h

        include/PhysiK/Math/AABB.h
        include/PhysiK/Math/Mat3.h
        include/PhysiK/Math/Quaternion.h
        include/PhysiK/Math/SparseBlockMatrix.h
        include/PhysiK/Math/Transform.h
        include/PhysiK/Math/Vec3.h
        include/PhysiK/Math/Vec4.h

        include/PhysiK/PhysicsData/Material.h
        include/PhysiK/PhysicsData/Node.h
        include/PhysiK/PhysicsData/Contact.h
        include/PhysiK/PhysicsData/Tet.h

        include/PhysiK/Components/Component.h
        include/PhysiK/Components/ComponentExecutionPriority.h
        include/PhysiK/Components/CollisionComponent.h
        include/PhysiK/Components/CollisionSphereComponent.h
        include/PhysiK/Components/ScriptComponent.h
        include/PhysiK/Components/SurfaceExtractionComponent.h
        include/PhysiK/Components/SurfaceVisualComponent.h
        include/PhysiK/Components/TetMeshComponent.h
        include/PhysiK/Components/TetMeshPhysicsComponent.h
        include/PhysiK/Components/TetMeshMapperComponent.h
        include/PhysiK/Components/TopologyMeshComponent.h
        include/PhysiK/Components/VisualMeshComponent.h
        include/PhysiK/Geometry/TetMeshGenerator.h

        include/PhysiK/Core/Collision/CollisionDetectionEngine.h
        include/PhysiK/Core/Events/EventSystem.h
        include/PhysiK/Core/Events/PhysicsEvent.h
        include/PhysiK/Core/Physics/PhysicsModel.h
        include/PhysiK/Core/Physics/FEM/FEMModel.h
        include/PhysiK/Core/PhysicsConnections/PhysicsConnection.h
        include/PhysiK/Core/PhysicsConnections/PointConnection.h
        include/PhysiK/Core/PhysicsConnections/SurfaceConnection.h
        include/PhysiK/Core/PhysicsConnections/LineConnection.h
        include/PhysiK/Core/PhysicsConnections/RigidBodyConnection.h
        include/PhysiK/Core/PhysicsConnections/RigidBodyOrientationConnection.h
        include/PhysiK/Core/Solvers/Linear/ConjugateGradientSolver.h
        include/PhysiK/Core/Solvers/Linear/DenseLinearSolver.h
        include/PhysiK/Core/Solvers/Linear/LinearSolver.h
        include/PhysiK/Core/Solvers/SolverData.h
        include/PhysiK/Core/World/World.h
    PRIVATE
        src/API/PhysiKAPI.cpp
        src/Components/CollisionComponent.cpp
        src/Components/CollisionSphereComponent.cpp
        src/Components/ScriptComponent.cpp
        src/Components/SurfaceExtractionComponent.cpp
        src/Components/SurfaceVisualComponent.cpp
        src/Components/TetMeshComponent.cpp
        src/Components/TetMeshPhysicsComponent.cpp
        src/Components/TetMeshMapperComponent.cpp
        src/Components/TopologyMeshComponent.cpp
        src/Components/VisualMeshComponent.cpp
        src/Geometry/TetMeshGenerator.cpp
        src/Core/Collision/CollisionDetectionEngine.cpp
        src/Core/Events/EventSystem.cpp
        src/Core/Physics/FEM/FEMModel.cpp
        src/Core/PhysicsConnections/PointConnection.cpp
        src/Core/PhysicsConnections/SurfaceConnection.cpp
        src/Core/PhysicsConnections/LineConnection.cpp
        src/Core/PhysicsConnections/RigidBodyConnection.cpp
        src/Core/PhysicsConnections/RigidBodyOrientationConnection.cpp
        src/Core/Solvers/Linear/ConjugateGradientSolver.cpp
        src/Core/Solvers/Linear/DenseLinearSolver.cpp
        src/Core/Solvers/Linear/LinearSolver.cpp
        src/Core/Solvers/SolverData.cpp
        src/Math/SparseBlockMatrix.cpp
        src/Core/World/World.cpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}/include" PREFIX "Public Headers" FILES
    include/PhysiK/API/Handles.h
    include/PhysiK/API/PhysiKAPI.h
    include/PhysiK/Math/AABB.h
    include/PhysiK/Math/Mat3.h
    include/PhysiK/Math/Quaternion.h
    include/PhysiK/Math/SparseBlockMatrix.h
    include/PhysiK/Math/Transform.h
    include/PhysiK/Math/Vec3.h
    include/PhysiK/Math/Vec4.h
    include/PhysiK/PhysicsData/Material.h
    include/PhysiK/PhysicsData/Node.h
    include/PhysiK/PhysicsData/Contact.h
    include/PhysiK/PhysicsData/Tet.h
    include/PhysiK/Components/Component.h
    include/PhysiK/Components/ComponentExecutionPriority.h
    include/PhysiK/Components/CollisionComponent.h
    include/PhysiK/Components/CollisionSphereComponent.h
    include/PhysiK/Components/ScriptComponent.h
    include/PhysiK/Components/SurfaceExtractionComponent.h
    include/PhysiK/Components/SurfaceVisualComponent.h
    include/PhysiK/Components/TetMeshComponent.h
    include/PhysiK/Components/TetMeshPhysicsComponent.h
    include/PhysiK/Components/TetMeshMapperComponent.h
    include/PhysiK/Components/TopologyMeshComponent.h
    include/PhysiK/Components/VisualMeshComponent.h
    include/PhysiK/Geometry/TetMeshGenerator.h
    include/PhysiK/Core/Collision/CollisionDetectionEngine.h
    include/PhysiK/Core/Events/EventSystem.h
    include/PhysiK/Core/Events/PhysicsEvent.h
    include/PhysiK/Core/Performance/PerformanceLogger.h
    include/PhysiK/Core/Physics/PhysicsModel.h
    include/PhysiK/Core/Physics/FEM/FEMModel.h
    include/PhysiK/Core/PhysicsConnections/PhysicsConnection.h
    include/PhysiK/Core/PhysicsConnections/PointConnection.h
    include/PhysiK/Core/PhysicsConnections/SurfaceConnection.h
    include/PhysiK/Core/PhysicsConnections/LineConnection.h
    include/PhysiK/Core/PhysicsConnections/RigidBodyConnection.h
    include/PhysiK/Core/PhysicsConnections/RigidBodyOrientationConnection.h
    include/PhysiK/Core/Solvers/Linear/ConjugateGradientSolver.h
    include/PhysiK/Core/Solvers/Linear/DenseLinearSolver.h
    include/PhysiK/Core/Solvers/Linear/LinearSolver.h
    include/PhysiK/Core/Solvers/SolverData.h
    include/PhysiK/Core/World/World.h
)

target_compile_definitions(PhysiK PRIVATE PHYSIK_BUILDING_DLL)

if(PHYSIK_ENABLE_PERF_LOGGING)
    target_compile_definitions(PhysiK PUBLIC PHYSIK_ENABLE_PERF_LOGGING)
    target_sources(
        PhysiK
        PUBLIC
            include/PhysiK/Core/Performance/PerformanceLogger.h
        PRIVATE
            src/Core/Performance/PerformanceLogger.cpp
    )
endif()

if(PHYSIK_BUILD_TESTS OR PHYSIK_BUILD_BENCHMARKS)
    target_compile_definitions(PhysiK PUBLIC PHYSIK_ENABLE_SOLVER_PROFILING)
endif()

if(PHYSIK_BUILD_BENCHMARKS)
    add_executable(PhysiKStressTest benchmarks/PhysiKStressTest.cpp)
    target_link_libraries(PhysiKStressTest PRIVATE PhysiK::PhysiK)
    target_compile_features(PhysiKStressTest PRIVATE cxx_std_17)

    if(WIN32)
        add_custom_command(
            TARGET PhysiKStressTest
            POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_FILE:PhysiK>
                $<TARGET_FILE_DIR:PhysiKStressTest>
        )
    endif()

endif()

if(PHYSIK_BUILD_TESTS)
    enable_testing()
    add_executable(PhysiKTests tests/PointConnectionTests.cpp)
    target_link_libraries(PhysiKTests PRIVATE PhysiK::PhysiK)
    target_compile_features(PhysiKTests PRIVATE cxx_std_17)

    add_executable(
        ComponentExecutionPriorityTests
        tests/ComponentExecutionPriorityTests.cpp
        src/Components/ScriptComponent.cpp
        src/Core/PhysicsConnections/PointConnection.cpp
        src/Core/World/World.cpp)
    target_link_libraries(ComponentExecutionPriorityTests PRIVATE PhysiK::PhysiK)
    target_compile_features(ComponentExecutionPriorityTests PRIVATE cxx_std_17)

    add_executable(
        ScriptComponentTests
        tests/ScriptComponentTests.cpp
        src/Components/ScriptComponent.cpp)
    target_link_libraries(ScriptComponentTests PRIVATE PhysiK::PhysiK)
    target_compile_features(ScriptComponentTests PRIVATE cxx_std_17)

    if(WIN32)
        add_custom_command(
            TARGET PhysiKTests
            POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_FILE:PhysiK>
                $<TARGET_FILE_DIR:PhysiKTests>
        )

        add_custom_command(
            TARGET ComponentExecutionPriorityTests
            POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_FILE:PhysiK>
                $<TARGET_FILE_DIR:ComponentExecutionPriorityTests>
        )

        add_custom_command(
            TARGET ScriptComponentTests
            POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_FILE:PhysiK>
                $<TARGET_FILE_DIR:ScriptComponentTests>
        )
    endif()

    add_test(NAME PhysiKTests COMMAND PhysiKTests)
    add_test(NAME ComponentExecutionPriorityTests COMMAND ComponentExecutionPriorityTests)
    add_test(NAME ScriptComponentTests COMMAND ScriptComponentTests)
endif()
