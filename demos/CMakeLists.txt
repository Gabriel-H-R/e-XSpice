# ============================================================
# Demo circuits — end-to-end smoke tests through the CLI.
#
# Each entry is "<netlist>|<signature line>". The signature pins the
# analysis that ran; analyses with adaptive stepping cannot pin an exact
# CSV row count, so the console line is the robust cross-platform check.
# ============================================================

set(DEMO_CASES
    "rc_step|Transient analysis completed"
    "ac_rlc_resonance|AC analysis completed"
    "noise_divider|Noise analysis completed"
    "sens_divider|Sensitivity analysis completed"
    "hb_clipper|Harmonic balance completed"
    "clipper_tran|Transient analysis completed"
    "bsim4_inverter|Transient analysis completed"
)

foreach(case IN LISTS DEMO_CASES)
    string(REPLACE "|" ";" _parts "${case}")
    list(GET _parts 0 _name)
    list(GET _parts 1 _pattern)
    add_test(NAME Demo_${_name}
             COMMAND ${CMAKE_COMMAND}
                     -DBINARY=$<TARGET_FILE:e-XSpice>
                     -DNETLIST=${CMAKE_CURRENT_SOURCE_DIR}/${_name}.cir
                     -DOUT=${CMAKE_CURRENT_BINARY_DIR}/${_name}.csv
                     -DPATTERN=${_pattern}
                     -P ${CMAKE_CURRENT_SOURCE_DIR}/check_demo.cmake)
endforeach()
