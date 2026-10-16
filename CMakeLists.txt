cmake_minimum_required(VERSION 3.15)

project(e-XSpice
    VERSION 1.0.0
    LANGUAGES CXX
)

# ============================================================
# C++ STANDARD
# ============================================================

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


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
# BSIM4 PHYSICS KERNEL
# ============================================================
# reference (default): pure-C++ subset, always builds, no extra headers.
# berkeley: upstream src/models/bsim4/code (b4ld.c/b4temp.c) via the
# adapter skeleton in src/models/external/. Requires SPICE3 shim headers
# (spice.h, cktdefs.h, devdefs.h — absent), so OFF by default.

option(EXSPICE_USE_BERKELEY_BSIM4 "Compile Berkeley BSIM4 C kernel (requires shim headers)" OFF)

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
message(STATUS "C++ Standard : C++${CMAKE_CXX_STANDARD}")
message(STATUS "Build Type   : ${CMAKE_BUILD_TYPE}")
message(STATUS "MPI          : ${MPI_FOUND}")
message(STATUS "OpenMP       : ${OpenMP_FOUND}")
message(STATUS "BSIM4 kernel : ${EXSPICE_USE_BERKELEY_BSIM4} (ON=Berkeley C, OFF=reference C++)")
message(STATUS "============================================")
message(STATUS "")