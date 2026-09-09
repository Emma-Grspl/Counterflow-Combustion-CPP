cmake_minimum_required(VERSION 3.20)

project(
    CounterflowCombustion
    VERSION 0.1.0
    DESCRIPTION "C++ solver for methane counterflow combustion"
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 3.4 REQUIRED NO_MODULE)

# ============================================================
# Core numerical library
# ============================================================

add_library(
    counterflow_core
    src/Grid.cpp
    src/BoundaryConditions.cpp
    src/NavierStokes.cpp
    src/PressurePoisson.cpp
)

target_include_directories(
    counterflow_core
    PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(
    counterflow_core
    PUBLIC
    Eigen3::Eigen
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(
        counterflow_core
        PRIVATE
        -Wall
        -Wextra
        -Wpedantic
    )
endif()

# ============================================================
# Main executable
# ============================================================

add_executable(
    counterflow_combustion
    src/main.cpp
)

target_link_libraries(
    counterflow_combustion
    PRIVATE
    counterflow_core
)

# ============================================================
# Tests
# ============================================================

enable_testing()

add_executable(
    test_config
    tests/test_config.cpp
)

target_link_libraries(
    test_config
    PRIVATE
    counterflow_core
)

add_test(
    NAME config
    COMMAND test_config
)

add_executable(
    test_grid
    tests/test_grid.cpp
)

target_link_libraries(
    test_grid
    PRIVATE
    counterflow_core
)

add_test(
    NAME grid
    COMMAND test_grid
)

add_executable(
    test_boundary_conditions
    tests/test_boundary_conditions.cpp
)

target_link_libraries(
    test_boundary_conditions
    PRIVATE
    counterflow_core
)

add_test(
    NAME boundary_conditions
    COMMAND test_boundary_conditions
)

add_executable(
    test_navier_stokes
    tests/test_navier_stokes.cpp
)

target_link_libraries(
    test_navier_stokes
    PRIVATE
    counterflow_core
)

add_test(
    NAME navier_stokes
    COMMAND test_navier_stokes
)

add_executable(
    test_pressure_poisson
    tests/test_pressure_poisson.cpp
)

target_link_libraries(
    test_pressure_poisson
    PRIVATE
    counterflow_core
)

add_test(
    NAME pressure_poisson
    COMMAND test_pressure_poisson
)
