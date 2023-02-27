# This is modified from the source of the GTest CMake module shipped with CMake but instead operates on
# UnitTest++ test drivers which expose --list-tests and --test-filter options.
# The GTest Cmake module is licensed under BSD-3 Clause

cmake_minimum_required(VERSION ${CMAKE_VERSION})

set(flush_tests_MODE WRITE)

macro(flush_script)
    file(${flush_tests_MODE} "${_CTEST_FILE}" "${script}")
    set(flush_tests_MODE APPEND)
    set(script "")
endmacro()

# Flushes tests_buffer to tests
macro(flush_tests_buffer)
    list(APPEND tests "${tests_buffer}")
    set(tests_buffer "")
endmacro()

macro(add_command NAME)
    set(_args "")
    foreach (_arg ${ARGN})
        if (_arg MATCHES "[^-./:a-zA-Z0-9_]")
            string(APPEND _args " [==[${_arg}]==]")
        else ()
            string(APPEND _args " ${_arg}")
        endif ()
    endforeach ()
    string(APPEND script "${NAME}(${_args})\n")
    string(LENGTH "${script}" _script_len)
    if (${_script_len} GREATER "50000")
        flush_script()
    endif ()
    # Unsets macro local variables to prevent leakage outside of this macro.
    unset(_args)
    unset(_script_len)
endmacro()

function(UnitTestPP_discover_tests_impl)

    cmake_parse_arguments(
            ""
            ""
            "TEST_EXECUTABLE;TEST_WORKING_DIR;TEST_PREFIX;TEST_SUFFIX;TEST_LIST;CTEST_FILE;TEST_DISCOVERY_TIMEOUT;TEST_XML_OUTPUT_DIR;TEST_FILTER"
            "TEST_EXTRA_ARGS;TEST_PROPERTIES;TEST_EXECUTOR"
            ${ARGN}
    )

    set(prefix "${_TEST_PREFIX}")
    set(suffix "${_TEST_SUFFIX}")
    set(extra_args "${_TEST_EXTRA_ARGS}")
    set(properties "${_TEST_PROPERTIES}")
    set(script)
    set(suite)
    set(tests)
    set(tests_buffer)

    if(_TEST_FILTER)
        set(filter "--test-filter=${_TEST_FILTER}")
    else()
        set(filter)
    endif()

    execute_process(
            COMMAND ${_TEST_EXECUTOR} "${_TEST_EXECUTABLE}" --list-tests ${filter}
            WORKING_DIRECTORY "${_TEST_WORKING_DIR}"
            TIMEOUT ${_TEST_DISCOVERY_TIMEOUT}
            OUTPUT_VARIABLE output
            RESULT_VARIABLE result
    )

    if(NOT ${result} EQUAL 0)
        message(FATAL_ERROR "error running test executable ${_TEST_EXECUTABLE}")
    endif()

    string(REPLACE [[;]] [[\;]] output "${output}")
    string(REPLACE "\n" ";" output "${output}")

    foreach(line ${output})
        set(testname "${line}")
        string(REPLACE [[\]] [[\\]] testname "${testname}")
        string(REPLACE [[;]] [[\;]] testname "${testname}")
        string(REPLACE [[$]] [[\$]] testname "${testname}")

        add_command(add_test
                "${testname}"
                ${_TEST_EXECUTOR}
                "${_TEST_EXECUTABLE}"
                "--test-filter=${testname}"
                ${extra_args}
                )
        add_command(set_tests_properties
                "${testname}"
                PROPERTIES
                    WORKING_DIRECTORY "${_TEST_WORKING_DIR}"
                    SKIP_REGULAR_EXPRESSION "\\\\[  SKIPPED\\\\]"
                    ${properties}
                )
        list(APPEND tests_buffer "${testname}")
        list(LENGTH tests_buffer tests_buffer_length)
        if (${tests_buffer_length} GREATER 250)
            flush_tests_buffer()
        endif()
    endforeach()

    flush_tests_buffer()
    add_command(set ${_TEST_LIST} ${tests})

    flush_script()

endfunction()


if(CMAKE_SCRIPT_MODE_FILE)
    UnitTestPP_discover_tests_impl(
            TEST_EXECUTABLE ${TEST_EXECUTABLE}
            TEST_EXECUTOR ${TEST_EXECUTOR}
            TEST_WORKING_DIR ${TEST_WORKING_DIR}
            TEST_PREFIX ${TEST_PREFIX}
            TEST_SUFFIX ${TEST_SUFFIX}
            TEST_FILTER ${TEST_FILTER}
            TEST_LIST ${TEST_LIST}
            CTEST_FILE ${CTEST_FILE}
            TEST_DISCOVERY_TIMEOUT ${TEST_DISCOVERY_TIMEOUT}
            TEST_XML_OUTPUT_DIR ${TEST_XML_OUTPUT_DIR}
            TEST_EXTRA_ARGS ${TEST_EXTRA_ARGS}
            TEST_PROPERTIES ${TEST_PROPERTIES}
   )

endif()
