# Tests CMakeLists.txt

# Simple test example
add_executable(test_sparse_matrix test_sparse_matrix.cpp)
target_link_libraries(test_sparse_matrix circuitsim_lib)

add_executable(test_resistor test_resistor.cpp)
target_link_libraries(test_resistor circuitsim_lib)

add_executable(test_dc_analysis test_dc_analysis.cpp)
target_link_libraries(test_dc_analysis circuitsim_lib)

# Register tests with CTest
add_test(NAME SparseMatrix COMMAND test_sparse_matrix)
add_test(NAME Resistor COMMAND test_resistor)
add_test(NAME DCAnalysis COMMAND mpirun -np 1 test_dc_analysis)
