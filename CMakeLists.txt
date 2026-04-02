cmake_minimum_required(VERSION 3.15)
project(CircuitSim VERSION 1.0.0 LANGUAGES CXX)

# C++17 is required for modern features
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Optimization and debug flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra")

# Find MPI (required for distributed parallelization)
find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_DIRS})

# Find OpenMP (for shared-memory parallelization)
find_package(OpenMP REQUIRED)

# Include project directories
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src
)

# Collect all source files
file(GLOB_RECURSE CORE_SOURCES 
    "src/core/*.cpp"
)

file(GLOB_RECURSE DEVICE_SOURCES 
    "src/devices/*.cpp"
)

file(GLOB_RECURSE ANALYSIS_SOURCES 
    "src/analysis/*.cpp"
)

file(GLOB_RECURSE SOLVER_SOURCES 
    "src/solver/*.cpp"
)

file(GLOB_RECURSE PARALLEL_SOURCES 
    "src/parallel/*.cpp"
)

# Create static library with all components
add_library(circuitsim_lib STATIC
    ${CORE_SOURCES}
    ${DEVICE_SOURCES}
    ${ANALYSIS_SOURCES}
    ${SOLVER_SOURCES}
    ${PARALLEL_SOURCES}
)

# Link MPI and OpenMP to the library
target_link_libraries(circuitsim_lib 
    PUBLIC 
    ${MPI_CXX_LIBRARIES}
    OpenMP::OpenMP_CXX
)

target_compile_options(circuitsim_lib 
    PUBLIC 
    ${MPI_CXX_COMPILE_FLAGS}
)

# Main executable
add_executable(circuitsim src/main.cpp)
target_link_libraries(circuitsim circuitsim_lib)

# Option to compile tests
option(BUILD_TESTS "Build test programs" ON)
if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

# Status messages
message(STATUS "MPI C++ Compiler: ${MPI_CXX_COMPILER}")
message(STATUS "OpenMP Flags: ${OpenMP_CXX_FLAGS}")
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
