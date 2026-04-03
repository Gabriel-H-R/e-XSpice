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

# Registrando os testes no CTest
add_test(NAME Divider_Test COMMAND test_divider)
add_test(NAME MultiLoop_Test COMMAND test_multiloop)
add_test(NAME RC_Step_Test COMMAND test_rc_step)
add_test(NAME RLC_Series_Test COMMAND test_rlc_series)