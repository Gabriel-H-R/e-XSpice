# check_demo.cmake — run one showcase netlist through the CLI and pin a
# signature line of its console output.
#
#   cmake -DBINARY=<e-XSpice> -DNETLIST=<file.cir> [-DOUT=<file.csv>]
#         -DPATTERN=<regex> -P check_demo.cmake

if(NOT DEFINED BINARY OR NOT DEFINED NETLIST OR NOT DEFINED PATTERN)
    message(FATAL_ERROR "check_demo.cmake: BINARY, NETLIST and PATTERN are required")
endif()

set(_cmd "${BINARY}" file "${NETLIST}")
if(DEFINED OUT AND OUT)
    list(APPEND _cmd "${OUT}")
endif()

execute_process(
    COMMAND ${_cmd}
    RESULT_VARIABLE _rc
    OUTPUT_VARIABLE _out
    ERROR_VARIABLE _err)

if(NOT _rc EQUAL 0)
    message(FATAL_ERROR "demo failed on ${NETLIST} (rc=${_rc})\n${_err}")
endif()

string(REGEX MATCH "${PATTERN}" _m "${_out}${_err}")
if(NOT _m)
    message(FATAL_ERROR "demo output does not match '${PATTERN}':\n${_out}\n${_err}")
endif()

message(STATUS "demo ok: ${NETLIST}")
