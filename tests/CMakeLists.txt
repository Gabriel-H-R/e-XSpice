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

# 8b. BJT level-2 (Gummel-Poon: qb, ISE/ISC, RB/RC/RE, substrate)
add_executable(test_bjt2 test_bjt2.cpp)
target_link_libraries(test_bjt2 e-XSpice_lib)

# 8c. MOSFET NMOS/PMOS level-1 (Shichman-Hodges + Meyer)
add_executable(test_mosfet test_mosfet.cpp)
target_link_libraries(test_mosfet e-XSpice_lib)

# 8d. MOSFET BSIM4 level-2 (NMOS2/PMOS2 cards, wrapper integration + MPI)
add_executable(test_bsim4 test_bsim4.cpp)
target_link_libraries(test_bsim4 e-XSpice_lib)

# 8e. DC sweep (.DC LIN/DEC/OCT over a V/I source, warm-started chain)
add_executable(test_dc_sweep test_dc_sweep.cpp)
target_link_libraries(test_dc_sweep e-XSpice_lib)

# 8f. Sparse complex LU (the AC small-signal factorization kernel)
add_executable(test_complex_lu test_complex_lu.cpp)
target_link_libraries(test_complex_lu e-XSpice_lib)

# 8g. AC small-signal analysis (.AC)
add_executable(test_ac test_ac.cpp)
target_link_libraries(test_ac e-XSpice_lib)

# 8h. Noise analysis (.NOISE)
add_executable(test_noise test_noise.cpp)
target_link_libraries(test_noise e-XSpice_lib)

# 8i. Sensitivity analyses (.SENS: DC, AC, transient)
add_executable(test_sens test_sens.cpp)
target_link_libraries(test_sens e-XSpice_lib)

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
add_test(NAME BJT_GummelPoon COMMAND test_bjt2)
add_test(NAME MOSFET_Level1 COMMAND test_mosfet)
add_test(NAME BSIM4_MOSFET COMMAND test_bsim4)
add_test(NAME DC_Sweep COMMAND test_dc_sweep)
add_test(NAME Complex_LU COMMAND test_complex_lu)
add_test(NAME AC_Analysis COMMAND test_ac)
add_test(NAME Noise_Analysis COMMAND test_noise)
add_test(NAME Sensitivity_Analysis COMMAND test_sens)

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
        # BSIM4 level-2 under MPI: exercises Mosfet::beginSolve/storeState and
        # the limiter vote on the BSIM4 branch on every rank (former 0xC0000005).
        add_test(NAME BSIM4_MOSFET_4Ranks
                 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 $<TARGET_FILE:test_bsim4>)
    endif()
    # Default-threshold mode (no EXSPICE_DISTRIBUTE override): exercises the
    # automatic assembly/solve selection under MPI.
    add_test(NAME MPI_Ladder_Auto_2Ranks
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:test_mpi_ladder> auto)
    # DC sweep under MPI: same deterministic sweep on every rank.
    add_test(NAME DC_Sweep_2Ranks
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:test_dc_sweep>)
endif()

# ============================================================
# End-to-end .DC netlist runs through the CLI (parser -> analysis -> CSV).
# ============================================================

# Diode clipper swept -10 V .. 15 V in 1 V steps (Xyce Users' Guide 1.4):
# 26 operating points, node columns named after the netlist nodes.
add_test(NAME DC_Sweep_Netlist
         COMMAND ${CMAKE_COMMAND}
                 -DBINARY=$<TARGET_FILE:e-XSpice>
                 -DNETLIST=${CMAKE_CURRENT_SOURCE_DIR}/dc_sweep_test_netlist.cir
                 -DEXPECT_ROWS=26
                 -DEXPECT_HEADER=VIN,Vin,Vout
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/check_csv.cmake)

# Nested sweep on a BSIM4 CMOS inverter: VDD (outer, 5 points) x VG (inner,
# 5 points) = 25 operating points, two swept columns in the CSV.
add_test(NAME DC_Sweep_Netlist_Nested
         COMMAND ${CMAKE_COMMAND}
                 -DBINARY=$<TARGET_FILE:e-XSpice>
                 -DNETLIST=${CMAKE_CURRENT_SOURCE_DIR}/dc_sweep_nested_netlist.cir
                 -DEXPECT_ROWS=25
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/check_csv.cmake)

# Flat sweep on a BSIM4 CMOS inverter: VG 0 -> 3 V in 0.5 V steps
# (7 operating points, one swept column).
add_test(NAME DC_Sweep_Netlist_Flat
         COMMAND ${CMAKE_COMMAND}
                 -DBINARY=$<TARGET_FILE:e-XSpice>
                 -DNETLIST=${CMAKE_CURRENT_SOURCE_DIR}/dc_sweep_flat_netlist.cir
                 -DEXPECT_ROWS=7
                 -DEXPECT_HEADER=VG,Vvdd,Vvg,Vvout
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/check_csv.cmake)

# Temperature sweep (.DC TEMP -40 125 25 on a copper resistor model):
# 8 points, one swept column (TEMP).
add_test(NAME DC_Sweep_Netlist_Temp
         COMMAND ${CMAKE_COMMAND}
                 -DBINARY=$<TARGET_FILE:e-XSpice>
                 -DNETLIST=${CMAKE_CURRENT_SOURCE_DIR}/dc_sweep_temp_netlist.cir
                 -DEXPECT_ROWS=7
                 -DEXPECT_HEADER=TEMP,Vin,Vout
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/check_csv.cmake)

# Same netlists under MPI: the distributed CSV must be byte-identical to the
# serial one (deterministic sweep on every rank).
if(MPIEXEC_EXECUTABLE)
    add_test(NAME DC_Sweep_Netlist_MPI
             COMMAND ${CMAKE_COMMAND}
                     -DBINARY=$<TARGET_FILE:e-XSpice>
                     -DNETLIST=${CMAKE_CURRENT_SOURCE_DIR}/dc_sweep_test_netlist.cir
                     -DEXPECT_ROWS=26
                     -DEXPECT_HEADER=VIN,Vin,Vout
                     -DMPIEXEC=${MPIEXEC_EXECUTABLE}
                     -DMPI_NP_FLAG=${MPIEXEC_NUMPROC_FLAG}
                     -DMPI_NP=2
                     -P ${CMAKE_CURRENT_SOURCE_DIR}/check_csv.cmake)
    add_test(NAME DC_Sweep_Netlist_Nested_MPI
             COMMAND ${CMAKE_COMMAND}
                     -DBINARY=$<TARGET_FILE:e-XSpice>
                     -DNETLIST=${CMAKE_CURRENT_SOURCE_DIR}/dc_sweep_nested_netlist.cir
                     -DEXPECT_ROWS=25
                     -DMPIEXEC=${MPIEXEC_EXECUTABLE}
                     -DMPI_NP_FLAG=${MPIEXEC_NUMPROC_FLAG}
                     -DMPI_NP=2
                     -P ${CMAKE_CURRENT_SOURCE_DIR}/check_csv.cmake)
    # `.OPTIONS SWEEPPAR`: sweep POINTS split across ranks (each rank keeps
    # its own continuation chain). Smoke test that the distributed-point path
    # completes; numeric equivalence to the sequential run is checked inside
    # test_dc_sweep (section 10), where a tolerance comparison is possible.
    add_test(NAME DC_Sweep_Netlist_SweepPar
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
                     $<TARGET_FILE:e-XSpice> file
                     ${CMAKE_CURRENT_SOURCE_DIR}/dc_sweep_sweeppar_netlist.cir
                     ${CMAKE_CURRENT_BINARY_DIR}/dcsweep_sweeppar.csv)
    set_tests_properties(DC_Sweep_Netlist_SweepPar PROPERTIES
        PASS_REGULAR_EXPRESSION "points distributed over 2 ranks")
endif()

# ================================================================
# AC small-signal analysis (.AC) — end to end through the CLI
# ================================================================
# RC low-pass, DEC 100 points/decade from 10 Hz to 1 MHz = 501 frequencies.
add_test(NAME AC_Netlist
         COMMAND ${CMAKE_COMMAND}
                 -DBINARY=$<TARGET_FILE:e-XSpice>
                 -DNETLIST=${CMAKE_CURRENT_SOURCE_DIR}/ac_rc_netlist.cir
                 -DEXPECT_ROWS=501
                 -DEXPECT_HEADER=frequency,Vin_mag,Vin_phase,Vout_mag,Vout_phase,V1_imag_mag,V1_imag_phase
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/check_csv.cmake)

# Series RLC (MNA with an inductor): the pivot-safe ordering path through the
# CLI, the branch-current columns of the CSV, and a resonant peak in the log.
add_test(NAME AC_Netlist_MNA
         COMMAND ${CMAKE_COMMAND}
                 -DBINARY=$<TARGET_FILE:e-XSpice>
                 -DNETLIST=${CMAKE_CURRENT_SOURCE_DIR}/ac_rlc_netlist.cir
                 -DEXPECT_ROWS=401
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/check_csv.cmake)

if(MPIEXEC_EXECUTABLE)
    # The AC sweep is deterministic by construction: frequencies are
    # independent and every rank factorizes the same replicated matrix, so the
    # distributed CSV must be BYTE-IDENTICAL to the serial one (unlike the
    # DC sweep, whose continuation chain is path dependent).
    add_test(NAME AC_Netlist_MPI
             COMMAND ${CMAKE_COMMAND}
                     -DBINARY=$<TARGET_FILE:e-XSpice>
                     -DNETLIST=${CMAKE_CURRENT_SOURCE_DIR}/ac_rc_netlist.cir
                     -DEXPECT_ROWS=501
                     -DEXPECT_HEADER=frequency,Vin_mag,Vin_phase,Vout_mag,Vout_phase,V1_imag_mag,V1_imag_phase
                     -DMPIEXEC=${MPIEXEC_EXECUTABLE}
                     -DMPI_NP_FLAG=${MPIEXEC_NUMPROC_FLAG}
                     -DMPI_NP=2
                     -P ${CMAKE_CURRENT_SOURCE_DIR}/check_csv.cmake)

    # `.OPTIONS ACFREQPAR`: frequencies split across ranks. Same guarantee —
    # bit-identical results, only the work distribution changes.
    add_test(NAME AC_Netlist_FreqPar
             COMMAND ${CMAKE_COMMAND}
                     -DBINARY=$<TARGET_FILE:e-XSpice>
                     -DNETLIST=${CMAKE_CURRENT_SOURCE_DIR}/ac_freqpar_netlist.cir
                     -DEXPECT_ROWS=41
                     -DMPIEXEC=${MPIEXEC_EXECUTABLE}
                     -DMPI_NP_FLAG=${MPIEXEC_NUMPROC_FLAG}
                     -DMPI_NP=2
                     -P ${CMAKE_CURRENT_SOURCE_DIR}/check_csv.cmake)
endif()

# ================================================================
# Sensitivity analyses (.SENS: DC, AC, transient) — end to end
# ================================================================
# The DC sensitivity writes one row per parameter (the exact closed form
# of a resistive divider: -2.5e-4 and +2.5e-4).
add_test(NAME Sens_DC_Netlist
         COMMAND ${CMAKE_COMMAND}
                 -DBINARY=$<TARGET_FILE:e-XSpice>
                 -DNETLIST=${CMAKE_CURRENT_SOURCE_DIR}/sens_dc_netlist.cir
                 -DEXPECT_ROWS=2
                 -DEXPECT_HEADER=parameter,sensitivity
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/check_csv.cmake)

# The AC sensitivity writes one row per frequency: the objective and the
# magnitude/phase of every parameter's sensitivity.
add_test(NAME Sens_AC_Netlist
         COMMAND ${CMAKE_COMMAND}
                 -DBINARY=$<TARGET_FILE:e-XSpice>
                 -DNETLIST=${CMAKE_CURRENT_SOURCE_DIR}/sens_ac_netlist.cir
                 -DEXPECT_ROWS=31
                 -DEXPECT_HEADER=frequency,obj_mag,obj_phase_deg,R1:R_mag,R1:R_phase_deg,C1:C_mag,C1:C_phase_deg
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/check_csv.cmake)

# The transient sensitivity writes one row per accepted time step.
add_test(NAME Sens_TRAN_Netlist
         COMMAND ${CMAKE_COMMAND}
                 -DBINARY=$<TARGET_FILE:e-XSpice>
                 -DNETLIST=${CMAKE_CURRENT_SOURCE_DIR}/sens_tran_netlist.cir
                 -DEXPECT_HEADER=time,obj_value,R1:R,C1:C
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/check_csv.cmake)

if(MPIEXEC_EXECUTABLE)
    # The AC sensitivity sweep is deterministic (the frequencies are
    # independent and every rank factorizes the same replicated matrix), so
    # the distributed CSV must be byte-identical to the serial one.
    add_test(NAME Sens_AC_Netlist_MPI
             COMMAND ${CMAKE_COMMAND}
                     -DBINARY=$<TARGET_FILE:e-XSpice>
                     -DNETLIST=${CMAKE_CURRENT_SOURCE_DIR}/sens_ac_netlist.cir
                     -DEXPECT_ROWS=31
                     -DEXPECT_HEADER=frequency,obj_mag,obj_phase_deg,R1:R_mag,R1:R_phase_deg,C1:C_mag,C1:C_phase_deg
                     -DMPIEXEC=${MPIEXEC_EXECUTABLE}
                     -DMPI_NP_FLAG=${MPIEXEC_NUMPROC_FLAG}
                     -DMPI_NP=2
                     -P ${CMAKE_CURRENT_SOURCE_DIR}/check_csv.cmake)
endif()

# ================================================================
# NOISE analysis (.NOISE) — end to end through the CLI
# ================================================================
# Resistive divider: the header and the per-generator columns (the .NOISE line
# carries pts_per_summary=1, which is what turns the breakdown on).
add_test(NAME Noise_Netlist
         COMMAND ${CMAKE_COMMAND}
                 -DBINARY=$<TARGET_FILE:e-XSpice>
                 -DNETLIST=${CMAKE_CURRENT_SOURCE_DIR}/noise_divider_netlist.cir
                 -DEXPECT_ROWS=61
                 -DEXPECT_HEADER=frequency,onoise_psd_v2_per_hz,inoise_psd_v2_per_hz,gain_mag,gain_phase_deg,R1_thermal_psd_v2_per_hz,R2_thermal_psd_v2_per_hz
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/check_csv.cmake)

if(MPIEXEC_EXECUTABLE)
    # The noise sweep splits frequencies, but every rank runs the SAME adjoint
    # solve on the same replicated matrix and the integral is a sum over the
    # whole sweep, so the distributed CSV must be byte-identical to the serial
    # one (same guarantee as the AC sweep).
    add_test(NAME Noise_Netlist_MPI
             COMMAND ${CMAKE_COMMAND}
                     -DBINARY=$<TARGET_FILE:e-XSpice>
                     -DNETLIST=${CMAKE_CURRENT_SOURCE_DIR}/noise_divider_netlist.cir
                     -DEXPECT_ROWS=61
                     -DEXPECT_HEADER=frequency,onoise_psd_v2_per_hz,inoise_psd_v2_per_hz,gain_mag,gain_phase_deg,R1_thermal_psd_v2_per_hz,R2_thermal_psd_v2_per_hz
                     -DMPIEXEC=${MPIEXEC_EXECUTABLE}
                     -DMPI_NP_FLAG=${MPIEXEC_NUMPROC_FLAG}
                     -DMPI_NP=2
                     -P ${CMAKE_CURRENT_SOURCE_DIR}/check_csv.cmake)
endif()











