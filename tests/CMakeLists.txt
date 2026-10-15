# CMakeLists.txt of tests's folder

# 1. Divider Circuit Test
add_executable(test_divider test_divider.cpp)
target_link_libraries(test_divider e-XSpice_lib)

# 2. Multi-Loop Circuit Test
add_executable(test_multiloop test_multiloop.cpp)
target_link_libraries(test_multiloop e-XSpice_lib)

# 3. RC Step Response Test
add_executable(test_rc_step test_rc_step.cpp)
target_link_libraries(test_rc_step e-XSpice_lib)

# 4. RLC Series Circuit Test
add_executable(test_rlc_series test_rlc_series.cpp)
target_link_libraries(test_rlc_series e-XSpice_lib)

# 5. Debug/verification harness (analytical checks on the DAE devices)
add_executable(debug_divider debug_divider.cpp)
target_link_libraries(debug_divider e-XSpice_lib)

# 6. Physics verification of the DAE device layer (analytic references)
add_executable(test_dae_physics test_dae_physics.cpp)
target_link_libraries(test_dae_physics e-XSpice_lib)

# 7. DeviceModel package (cards, registry, sharing, temperature)
add_executable(test_models test_models.cpp)
target_link_libraries(test_models e-XSpice_lib)

# 8. BJT NPN/PNP (Ebers-Moll transport model, charge partition)
add_executable(test_bjt test_bjt.cpp)
target_link_libraries(test_bjt e-XSpice_lib)

# DC convergence diagnostics: current-source biased BJT (limiting walk)
add_executable(debug_bjt_dc debug_bjt_dc.cpp)
target_link_libraries(debug_bjt_dc e-XSpice_lib)

# Registrando os testes no CTest
add_test(NAME Divider_Test COMMAND test_divider)
add_test(NAME MultiLoop_Test COMMAND test_multiloop)
add_test(NAME RC_Step_Test COMMAND test_rc_step)
add_test(NAME RLC_Series_Test COMMAND test_rlc_series)
add_test(NAME Debug_Divider COMMAND debug_divider)
add_test(NAME DAE_Physics COMMAND test_dae_physics)
add_test(NAME Device_Models COMMAND test_models)
add_test(NAME BJT_EbersMoll COMMAND test_bjt)

# RL step debug
add_executable(debug_rl debug_rl.cpp)
target_link_libraries(debug_rl e-XSpice_lib)
add_test(NAME Debug_RL COMMAND debug_rl)


# RLC shape debug
add_executable(debug_rlc debug_rlc.cpp)
target_link_libraries(debug_rlc e-XSpice_lib)


# Rectifier/pulse-edge debug
add_executable(debug_rectifier debug_rectifier.cpp)
target_link_libraries(debug_rectifier e-XSpice_lib)


# Benchmark teardown isolation
add_executable(debug_benchmark debug_benchmark.cpp)
target_link_libraries(debug_benchmark e-XSpice_lib)
add_test(NAME Benchmark_GroundedAnode COMMAND debug_benchmark)

# MPI equivalence: serial vs distributed on an RC ladder.
# Runs serially by default; runs under mpiexec when available.
add_executable(test_mpi_ladder test_mpi_ladder.cpp)
target_link_libraries(test_mpi_ladder e-XSpice_lib)
add_test(NAME MPI_Ladder_Serial COMMAND test_mpi_ladder)
if(MPIEXEC_EXECUTABLE)
    add_test(NAME MPI_Ladder_2Ranks
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:test_mpi_ladder>)
    if(MPIEXEC_MAX_NUMPROCS GREATER 3 OR MPIEXEC_MAX_NUMPROCS EQUAL "")
        add_test(NAME MPI_Ladder_4Ranks
                 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:test_mpi_ladder>)
    endif()
    # Default-threshold mode (no EXSPICE_DISTRIBUTE override): exercises the
    # automatic assembly/solve selection under MPI.
    add_test(NAME MPI_Ladder_Auto_2Ranks
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:test_mpi_ladder> auto)
endif()







