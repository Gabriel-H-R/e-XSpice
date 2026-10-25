cmake_minimum_required(VERSION 3.15)

project(e-XSpice
    VERSION 1.0.0
    LANGUAGES CXX C
)

# ============================================================
# C++ / C STANDARD
# ============================================================

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)


# ============================================================
# BUILD TYPE
# ============================================================

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()


# ============================================================
# COMPILER FLAGS
# ============================================================

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(-O3 -march=native -DNDEBUG)
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-g -O0 -Wall -Wextra)
endif()


# ============================================================
# MPI
# ============================================================

find_package(MPI)

if(MPI_FOUND)
    message(STATUS "MPI found.")
    message(STATUS "MPI C++ Compiler: ${MPI_CXX_COMPILER}")
else()
    message(STATUS "MPI not found. Continuing without MPI.")
endif()


# ============================================================
# OPENMP
# ============================================================

find_package(OpenMP)

if(OpenMP_FOUND)
    message(STATUS "OpenMP found.")
    message(STATUS "OpenMP Flags: ${OpenMP_CXX_FLAGS}")
else()
    message(STATUS "OpenMP not found. Continuing without OpenMP.")
endif()


# ============================================================
# INCLUDE DIRECTORIES
# ============================================================

include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/models/core_models
    ${CMAKE_SOURCE_DIR}/src/models/external          # SPICE shim headers
    ${CMAKE_SOURCE_DIR}/src/models/bsim4/code        # Berkeley BSIM4 headers
    ${CMAKE_SOURCE_DIR}/libs/symbolic
)


# ============================================================
# SOURCE FILES
# ============================================================

file(GLOB_RECURSE CORE_SOURCES
    "${CMAKE_SOURCE_DIR}/src/core/*.cpp"
)

file(GLOB_RECURSE DEVICE_SOURCES
    "${CMAKE_SOURCE_DIR}/src/devices/*.cpp"
)

file(GLOB_RECURSE ANALYSIS_SOURCES
    "${CMAKE_SOURCE_DIR}/src/analysis/*.cpp"
)

file(GLOB_RECURSE SOLVER_SOURCES
    "${CMAKE_SOURCE_DIR}/src/solver/*.cpp"
)

file(GLOB_RECURSE PARALLEL_SOURCES
    "${CMAKE_SOURCE_DIR}/src/parallel/*.cpp"
)

file(GLOB_RECURSE MODEL_SOURCES
    "${CMAKE_SOURCE_DIR}/src/models/*.cpp"
)

file(GLOB_RECURSE SYMBOLIC_SOURCES
    "${CMAKE_SOURCE_DIR}/libs/symbolic/*.cpp"
)

# Berkeley BSIM4 C sources (compiled as C, verbatim, no changes).
# Only the load/temp functions are needed; setup/getic/destroy and the
# parameter tables (b4.c) require the full SPICE environment and are
# handled manually by BSim4Adapter.
set(BERKELEY_BSIM4_SOURCES
    "${CMAKE_SOURCE_DIR}/src/models/bsim4/code/b4ld.c"
    "${CMAKE_SOURCE_DIR}/src/models/bsim4/code/b4temp.c"
)


# ============================================================
# STATIC LIBRARY
# ============================================================

add_library(e-XSpice_lib STATIC

    ${CORE_SOURCES}
    ${DEVICE_SOURCES}
    ${ANALYSIS_SOURCES}
    ${SOLVER_SOURCES}
    ${PARALLEL_SOURCES}
    ${MODEL_SOURCES}
    ${SYMBOLIC_SOURCES}
    ${BERKELEY_BSIM4_SOURCES}

)


# ============================================================
# MPI
# ============================================================

if(MPI_FOUND)

    target_link_libraries(e-XSpice_lib
        PUBLIC
        MPI::MPI_CXX
    )

endif()


# ============================================================
# OPENMP
# ============================================================

if(OpenMP_FOUND)

    target_link_libraries(e-XSpice_lib
        PUBLIC
        OpenMP::OpenMP_CXX
    )

endif()


# ============================================================
# MAIN EXECUTABLE
# ============================================================

add_executable(e-XSpice
    src/main.cpp
)

target_link_libraries(e-XSpice
    PRIVATE
    e-XSpice_lib
)


# ============================================================
# BSIM4 PHYSICS KERNEL (Berkeley C, now always available)
# ============================================================
# The SPICE shim headers (src/models/external/SpiceHeaders.h) provide the
# minimal type declarations needed to compile the Berkeley BSIM4 C code
# unchanged. The BSim4Adapter wraps BSIM4load() into e-XSpice's DAE contract.

option(EXSPICE_USE_BERKELEY_BSIM4 "Use Berkeley BSIM4 C kernel via BSim4Adapter" ON)

if(EXSPICE_USE_BERKELEY_BSIM4)
    add_compile_definitions(EXSPICE_USE_BERKELEY_BSIM4)
endif()


# ============================================================
# TESTS
# ============================================================

option(BUILD_TESTS "Build test programs" ON)

if(BUILD_TESTS)

    enable_testing()

    add_subdirectory(tests)

endif()


# ============================================================
# STATUS
# ============================================================

message(STATUS "")
message(STATUS "========== e-XSpice Configuration ==========")
message(STATUS "C++ Standard  : C++${CMAKE_CXX_STANDARD}")
message(STATUS "C Standard    : C${CMAKE_C_STANDARD}")
message(STATUS "Build Type    : ${CMAKE_BUILD_TYPE}")
message(STATUS "MPI           : ${MPI_FOUND}")
message(STATUS "OpenMP        : ${OpenMP_FOUND}")
message(STATUS "BSIM4 kernel  : ${EXSPICE_USE_BERKELEY_BSIM4} (ON=Berkeley C)")
message(STATUS "============================================")
message(STATUS "")
