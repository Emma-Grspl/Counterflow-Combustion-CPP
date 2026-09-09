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

add_executable(
    counterflow_combustion
    src/main.cpp
)

target_include_directories(
    counterflow_combustion
    PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(
        counterflow_combustion
        PRIVATE
        -Wall
        -Wextra
        -Wpedantic
    )
endif()

enable_testing()

add_executable(
    test_config
    tests/test_config.cpp
)

target_include_directories(
    test_config
    PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(
        test_config
        PRIVATE
        -Wall
        -Wextra
        -Wpedantic
    )
endif()

add_test(
    NAME config
    COMMAND test_config
)
