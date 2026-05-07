cmake_minimum_required(VERSION 3.20)

project(PhysiK VERSION 0.1.0 LANGUAGES CXX)

option(PHYSIK_BUILD_TESTS "Build PhysiK tests" ON)

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
    PRIVATE
        src/API/PhysiKAPI.cpp
        src/Core/World/World.cpp
)

target_compile_definitions(PhysiK PRIVATE PHYSIK_BUILDING_DLL)

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
endif()
