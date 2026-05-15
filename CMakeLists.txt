cmake_minimum_required(VERSION 3.20)

project(PhysiK VERSION 0.1.0 LANGUAGES CXX)

option(PHYSIK_BUILD_TESTS "Build PhysiK tests" ON)
option(PHYSIK_ENABLE_MKL "Enable optional Intel MKL linear solver support" OFF)

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
        include/PhysiK/Math/CSRMatrix.h
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
        include/PhysiK/Components/CollisionComponent.h
        include/PhysiK/Components/CollisionSphereComponent.h
        include/PhysiK/Components/TetMeshComponent.h

        include/PhysiK/Core/Collision/CollisionDetectionEngine.h
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
        include/PhysiK/Core/Solvers/Linear/MKLLinearSolver.h
        include/PhysiK/Core/Solvers/SolverData.h
        include/PhysiK/Core/World/World.h
    PRIVATE
        src/API/PhysiKAPI.cpp
        src/Components/CollisionComponent.cpp
        src/Components/CollisionSphereComponent.cpp
        src/Components/TetMeshComponent.cpp
        src/Core/Collision/CollisionDetectionEngine.cpp
        src/Core/Physics/FEM/FEMModel.cpp
        src/Core/PhysicsConnections/PointConnection.cpp
        src/Core/PhysicsConnections/SurfaceConnection.cpp
        src/Core/PhysicsConnections/LineConnection.cpp
        src/Core/PhysicsConnections/RigidBodyConnection.cpp
        src/Core/PhysicsConnections/RigidBodyOrientationConnection.cpp
        src/Core/Solvers/Linear/ConjugateGradientSolver.cpp
        src/Core/Solvers/Linear/DenseLinearSolver.cpp
        src/Core/Solvers/Linear/LinearSolver.cpp
        src/Core/Solvers/Linear/MKLLinearSolver.cpp
        src/Core/Solvers/SolverData.cpp
        src/Math/CSRMatrix.cpp
        src/Math/SparseBlockMatrix.cpp
        src/Core/World/World.cpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}/include" PREFIX "Public Headers" FILES
    include/PhysiK/API/Handles.h
    include/PhysiK/API/PhysiKAPI.h
    include/PhysiK/Math/AABB.h
    include/PhysiK/Math/CSRMatrix.h
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
    include/PhysiK/Components/CollisionComponent.h
    include/PhysiK/Components/CollisionSphereComponent.h
    include/PhysiK/Components/TetMeshComponent.h
    include/PhysiK/Core/Collision/CollisionDetectionEngine.h
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
    include/PhysiK/Core/Solvers/Linear/MKLLinearSolver.h
    include/PhysiK/Core/Solvers/SolverData.h
    include/PhysiK/Core/World/World.h
)

target_compile_definitions(PhysiK PRIVATE PHYSIK_BUILDING_DLL)

if(PHYSIK_ENABLE_MKL)
    find_package(MKL CONFIG QUIET)

    if(NOT MKL_FOUND)
        message(FATAL_ERROR
            "PHYSIK_ENABLE_MKL is ON, but Intel MKL was not found. "
            "Install Intel oneAPI MKL and set MKL_DIR or CMAKE_PREFIX_PATH "
            "so CMake can locate MKLConfig.cmake.")
    endif()

    if(NOT TARGET MKL::MKL)
        message(FATAL_ERROR
            "PHYSIK_ENABLE_MKL is ON and MKL was found, but the expected "
            "CMake target MKL::MKL is unavailable.")
    endif()

    target_compile_definitions(PhysiK PUBLIC PHYSIK_ENABLE_MKL)
    target_link_libraries(PhysiK PRIVATE MKL::MKL)

    set(PHYSIK_MKL_RUNTIME_SEARCH_DIRS)
    if(DEFINED MKL_DIR)
        get_filename_component(PHYSIK_MKL_ROOT_FROM_DIR "${MKL_DIR}/../../.." ABSOLUTE)
        get_filename_component(PHYSIK_ONEAPI_ROOT_FROM_MKL "${PHYSIK_MKL_ROOT_FROM_DIR}/../.." ABSOLUTE)
        list(APPEND PHYSIK_MKL_RUNTIME_SEARCH_DIRS
            "${PHYSIK_MKL_ROOT_FROM_DIR}/redist/intel64"
            "${PHYSIK_MKL_ROOT_FROM_DIR}/bin"
            "${PHYSIK_ONEAPI_ROOT_FROM_MKL}/compiler/latest/redist/intel64_win/compiler"
            "${PHYSIK_ONEAPI_ROOT_FROM_MKL}/compiler/latest/bin"
        )
        file(GLOB PHYSIK_MKL_VERSIONED_RUNTIME_DIRS
            "${PHYSIK_ONEAPI_ROOT_FROM_MKL}/*/bin"
            "${PHYSIK_ONEAPI_ROOT_FROM_MKL}/mkl/*/bin"
            "${PHYSIK_ONEAPI_ROOT_FROM_MKL}/compiler/*/bin"
        )
        list(APPEND PHYSIK_MKL_RUNTIME_SEARCH_DIRS
            ${PHYSIK_MKL_VERSIONED_RUNTIME_DIRS}
        )
    endif()
    if(DEFINED ENV{MKLROOT})
        list(APPEND PHYSIK_MKL_RUNTIME_SEARCH_DIRS
            "$ENV{MKLROOT}/redist/intel64"
            "$ENV{MKLROOT}/bin"
        )
    endif()
    if(DEFINED ENV{CMPLR_ROOT})
        list(APPEND PHYSIK_MKL_RUNTIME_SEARCH_DIRS
            "$ENV{CMPLR_ROOT}/redist/intel64_win/compiler"
            "$ENV{CMPLR_ROOT}/bin"
        )
    endif()
    if(DEFINED ENV{ONEAPI_ROOT})
        list(APPEND PHYSIK_MKL_RUNTIME_SEARCH_DIRS
            "$ENV{ONEAPI_ROOT}/mkl/latest/redist/intel64"
            "$ENV{ONEAPI_ROOT}/compiler/latest/redist/intel64_win/compiler"
            "$ENV{ONEAPI_ROOT}/compiler/latest/bin"
        )
    endif()
    list(APPEND PHYSIK_MKL_RUNTIME_SEARCH_DIRS
        "C:/Program Files (x86)/Intel/oneAPI/mkl/latest/redist/intel64"
        "C:/Program Files (x86)/Intel/oneAPI/mkl/latest/bin"
        "C:/Program Files (x86)/Intel/oneAPI/compiler/latest/redist/intel64_win/compiler"
        "C:/Program Files (x86)/Intel/oneAPI/compiler/latest/bin"
        "C:/Program Files/Intel/oneAPI/mkl/latest/redist/intel64"
        "C:/Program Files/Intel/oneAPI/mkl/latest/bin"
        "C:/Program Files/Intel/oneAPI/compiler/latest/redist/intel64_win/compiler"
        "C:/Program Files/Intel/oneAPI/compiler/latest/bin"
    )
    file(GLOB PHYSIK_COMMON_ONEAPI_RUNTIME_DIRS
        "C:/Program Files (x86)/Intel/oneAPI/*/bin"
        "C:/Program Files (x86)/Intel/oneAPI/mkl/*/bin"
        "C:/Program Files (x86)/Intel/oneAPI/compiler/*/bin"
        "C:/Program Files/Intel/oneAPI/*/bin"
        "C:/Program Files/Intel/oneAPI/mkl/*/bin"
        "C:/Program Files/Intel/oneAPI/compiler/*/bin"
    )
    list(APPEND PHYSIK_MKL_RUNTIME_SEARCH_DIRS
        ${PHYSIK_COMMON_ONEAPI_RUNTIME_DIRS}
    )
    list(REMOVE_DUPLICATES PHYSIK_MKL_RUNTIME_SEARCH_DIRS)

    set(PHYSIK_MKL_RUNTIME_DLLS)
    foreach(PHYSIK_MKL_RUNTIME_DLL_NAME
        mkl_core.3.dll
        mkl_intel_thread.3.dll
        libiomp5md.dll
    )
        string(MAKE_C_IDENTIFIER
            "${PHYSIK_MKL_RUNTIME_DLL_NAME}"
            PHYSIK_MKL_RUNTIME_DLL_ID
        )
        find_file(
            PHYSIK_MKL_RUNTIME_DLL_${PHYSIK_MKL_RUNTIME_DLL_ID}
            NAMES ${PHYSIK_MKL_RUNTIME_DLL_NAME}
            PATHS ${PHYSIK_MKL_RUNTIME_SEARCH_DIRS}
            NO_DEFAULT_PATH
        )
        if(NOT PHYSIK_MKL_RUNTIME_DLL_${PHYSIK_MKL_RUNTIME_DLL_ID})
            find_file(
                PHYSIK_MKL_RUNTIME_DLL_${PHYSIK_MKL_RUNTIME_DLL_ID}
                NAMES ${PHYSIK_MKL_RUNTIME_DLL_NAME}
            )
        endif()
        if(NOT PHYSIK_MKL_RUNTIME_DLL_${PHYSIK_MKL_RUNTIME_DLL_ID})
            message(FATAL_ERROR
                "PHYSIK_ENABLE_MKL is ON, but required MKL runtime DLL "
                "${PHYSIK_MKL_RUNTIME_DLL_NAME} was not found. Set MKLROOT, "
                "ONEAPI_ROOT, or CMPLR_ROOT so CMake can locate the oneAPI "
                "runtime redistributables.")
        endif()
        list(APPEND PHYSIK_MKL_RUNTIME_DLLS
            "${PHYSIK_MKL_RUNTIME_DLL_${PHYSIK_MKL_RUNTIME_DLL_ID}}"
        )
    endforeach()
    foreach(PHYSIK_MKL_OPTIONAL_RUNTIME_DLL_NAME
        mkl_intel_lp64.3.dll
        mkl_rt.3.dll
        mkl_def.3.dll
        mkl_mc3.3.dll
        mkl_avx2.3.dll
        mkl_avx512.3.dll
        mkl_avx10.3.dll
    )
        string(MAKE_C_IDENTIFIER
            "${PHYSIK_MKL_OPTIONAL_RUNTIME_DLL_NAME}"
            PHYSIK_MKL_OPTIONAL_RUNTIME_DLL_ID
        )
        find_file(
            PHYSIK_MKL_RUNTIME_DLL_${PHYSIK_MKL_OPTIONAL_RUNTIME_DLL_ID}
            NAMES ${PHYSIK_MKL_OPTIONAL_RUNTIME_DLL_NAME}
            PATHS ${PHYSIK_MKL_RUNTIME_SEARCH_DIRS}
            NO_DEFAULT_PATH
        )
        if(NOT PHYSIK_MKL_RUNTIME_DLL_${PHYSIK_MKL_OPTIONAL_RUNTIME_DLL_ID})
            find_file(
                PHYSIK_MKL_RUNTIME_DLL_${PHYSIK_MKL_OPTIONAL_RUNTIME_DLL_ID}
                NAMES ${PHYSIK_MKL_OPTIONAL_RUNTIME_DLL_NAME}
            )
        endif()
        if(PHYSIK_MKL_RUNTIME_DLL_${PHYSIK_MKL_OPTIONAL_RUNTIME_DLL_ID})
            list(APPEND PHYSIK_MKL_RUNTIME_DLLS
                "${PHYSIK_MKL_RUNTIME_DLL_${PHYSIK_MKL_OPTIONAL_RUNTIME_DLL_ID}}"
            )
        else()
            message(WARNING
                "Optional MKL runtime DLL ${PHYSIK_MKL_OPTIONAL_RUNTIME_DLL_NAME} "
                "was not found and will not be copied.")
        endif()
    endforeach()

    function(physik_copy_runtime_dlls target)
        if(WIN32)
            add_custom_command(
                TARGET ${target}
                POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    ${PHYSIK_MKL_RUNTIME_DLLS}
                    $<TARGET_FILE_DIR:${target}>
            )
        endif()
    endfunction()
endif()

if(PHYSIK_BUILD_TESTS)
    enable_testing()
    add_executable(PhysiKTests tests/PointConnectionTests.cpp)
    target_link_libraries(PhysiKTests PRIVATE PhysiK::PhysiK)
    target_compile_features(PhysiKTests PRIVATE cxx_std_17)

    if(WIN32)
        add_custom_command(
            TARGET PhysiKTests
            POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                $<TARGET_FILE:PhysiK>
                $<TARGET_FILE_DIR:PhysiKTests>
        )
    endif()

    add_test(NAME PhysiKTests COMMAND PhysiKTests)

    if(PHYSIK_ENABLE_MKL)
        add_executable(PhysiKMKLTests tests/MKLLinearSolverTests.cpp)
        target_link_libraries(PhysiKMKLTests PRIVATE PhysiK::PhysiK)
        target_compile_features(PhysiKMKLTests PRIVATE cxx_std_17)

        if(WIN32)
            add_custom_command(
                TARGET PhysiKMKLTests
                POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    $<TARGET_FILE:PhysiK>
                    $<TARGET_FILE_DIR:PhysiKMKLTests>
            )
        endif()
        physik_copy_runtime_dlls(PhysiKMKLTests)

        add_test(NAME PhysiKMKLTests COMMAND PhysiKMKLTests)
    endif()
endif()
