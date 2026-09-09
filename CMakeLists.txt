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

add_library(
    counterflow_core
    src/Grid.cpp
)

target_include_directories(
    counterflow_core
    PUBLIC
    ${PROJECT_SOURCE_DIR}/include
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

add_executable(
    counterflow_combustion
    src/main.cpp
)

target_link_libraries(
    counterflow_combustion
    PRIVATE
    counterflow_core
)

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
