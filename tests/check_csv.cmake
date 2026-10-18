cmake_policy(SET CMP0007 NEW)   # list() must not silently drop empty elements

# ============================================================
# check_csv.cmake — end-to-end validation of a .DC netlist run
#
#   cmake -DBINARY=<e-XSpice> -DNETLIST=<file.cir> -DEXPECT_ROWS=<n>
#         [-DEXPECT_HEADER=<csv header line>]
#         [-DMPIEXEC=<mpiexec> -DMPI_NP_FLAG=-n -DMPI_NP=2]
#         [-P check_csv.cmake]
#
# Runs the simulator on the netlist, then checks that
#   - the run succeeded,
#   - the CSV header and the number of data rows match,
#   - when MPIEXEC is given, the distributed run produces a CSV
#     byte-identical to the serial one (deterministic sweep).
# ============================================================

if(NOT DEFINED BINARY OR NOT DEFINED NETLIST)
    message(FATAL_ERROR "check_csv.cmake: BINARY and NETLIST are required")
endif()

set(_serial "${CMAKE_CURRENT_BINARY_DIR}/dcsweep_check_serial.csv")
set(_mpi    "${CMAKE_CURRENT_BINARY_DIR}/dcsweep_check_mpi.csv")

execute_process(
    COMMAND "${BINARY}" file "${NETLIST}" "${_serial}"
    RESULT_VARIABLE _rc OUTPUT_VARIABLE _out ERROR_VARIABLE _err)
if(NOT _rc EQUAL 0)
    message(FATAL_ERROR "e-XSpice failed on ${NETLIST} (rc=${_rc})\n${_err}")
endif()
if(NOT EXISTS "${_serial}")
    message(FATAL_ERROR "CSV was not produced: ${_serial}")
endif()

file(READ "${_serial}" _content)
string(REPLACE "\n" ";" _lines "${_content}")
list(LENGTH _lines _n)
# Trailing newline produces one empty trailing element.
math(EXPR _rows "${_n} - 1")
if(_rows LESS 1)
    message(FATAL_ERROR "CSV has no data rows: ${_serial}")
endif()

if(DEFINED EXPECT_HEADER)
    string(REGEX MATCH "^([^\n]*)\n" _m "${_content}")
    if(NOT CMAKE_MATCH_1 STREQUAL EXPECT_HEADER)
        message(FATAL_ERROR "CSV header mismatch:\n  got:      ${CMAKE_MATCH_1}\n  expected: ${EXPECT_HEADER}")
    endif()
endif()

if(DEFINED EXPECT_ROWS)
    # EXPECT_ROWS counts DATA rows (the header is not one of them).
    math(EXPR _data "${_rows} - 1")
    if(NOT _data EQUAL EXPECT_ROWS)
        message(FATAL_ERROR "CSV row count mismatch: got ${_data}, expected ${EXPECT_ROWS}")
    endif()
endif()

message(STATUS "CSV ok: ${NETLIST} (${_rows} lines, header '${CMAKE_MATCH_1}')")

if(DEFINED MPIEXEC AND MPIEXEC)
    execute_process(
        COMMAND "${MPIEXEC}" "${MPI_NP_FLAG}" "${MPI_NP}" "${BINARY}" file "${NETLIST}" "${_mpi}"
        RESULT_VARIABLE _rc2 OUTPUT_VARIABLE _out2 ERROR_VARIABLE _err2)
    if(NOT _rc2 EQUAL 0)
        message(FATAL_ERROR "mpiexec run failed (rc=${_rc2})\n${_err2}")
    endif()
    file(READ "${_mpi}" _content2)
    if(NOT _content STREQUAL _content2)
        message(FATAL_ERROR "Distributed CSV differs from the serial one: ${_mpi}")
    endif()
    message(STATUS "Serial/distributed CSV identical (${MPI_NP} ranks)")
endif()