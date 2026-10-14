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

# Registrando os testes no CTest
add_test(NAME Divider_Test COMMAND test_divider)
add_test(NAME MultiLoop_Test COMMAND test_multiloop)
add_test(NAME RC_Step_Test COMMAND test_rc_step)
add_test(NAME RLC_Series_Test COMMAND test_rlc_series)
add_test(NAME Debug_Divider COMMAND debug_divider)
add_test(NAME DAE_Physics COMMAND test_dae_physics)

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
