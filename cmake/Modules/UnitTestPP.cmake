# This is modified from the source of the GTest CMake module shipped with CMake but instead operates on
# UnitTest++ test drivers which expose --list-tests and --test-filter options.
# The GTest Cmake module is licensed under BSD-3 Clause

cmake_policy(PUSH)
cmake_policy(SET CMP0057 NEW)

function(UnitTestPP_discover_tests TARGET)
    cmake_parse_arguments(
            ""
            ""
            "TEST_PREFIX;TEST_SUFFIX;WORKING_DIRECTORY;TEST_LIST;DISCOVER_TIMEOUT;XML_OUTPUT_DIR;DISCOVERY_MODE"
            "EXTRA_ARGS;PROPERTIES;TEST_FILTER"
            ${ARGN}
            )

    if (NOT _WORKING_DIRECTORY)
        set(_WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
    endif()
    if (NOT _TEST_LIST)
        set(_TEST_LIST ${TARGET}_TESTS)
    endif()
    if(NOT _DISCOVERY_TIMEOUT)
        set(_DISCOVERY_TIMEOUT 5)
    endif()
    if (NOT _DISCOVERY_MODE)
        set(_DISCOVERY_MODE "POST_BUILD")
    endif()

    get_property(
            has_counter
            TARGET ${TARGET}
            PROPERTY CTEST_DISCOVERED_TEST_COUNTER
            SET
    )
    if (has_counter)
        get_property(
                counter
                TARGET ${TARGET}
                PROPERTY CTEST_DISCOVERED_TEST_COUNTER
        )
        math(EXPR counter "${counter} + 1")
    else()
        set(counter 1)
    endif()
    set_property(
            TARGET ${TARGET}
            PROPERTY CTEST_DISCOVERED_TEST_COUNTER
            ${counter}
    )

    set(ctest_file_base "${CMAKE_CURRENT_BINARY_DIR}/${TARGET}[${counter}]")
    set(ctest_include_file "${ctest_file_base}_include.cmake")
    set(ctest_tests_file "${ctest_file_base}_tests.cmake")

    get_property(crosscompiling_emulator
            TARGET ${TARGET}
            PROPERTY CROSSCOMPILING_EMULATOR
            )

    if (_DISCOVERY_MODE STREQUAL "POST_BUILD")
        add_custom_command(
                TARGET ${TARGET} POST_BUILD
                BYPRODUCTS "${ctest_tests_file}"
                COMMAND "${CMAKE_COMMAND}"
                        -D "TEST_TARGET=${TARGET}"
                        -D "TEST_EXECUTABLE=$<TARGET_FILE:${TARGET}>"
                        -D "TEST_EXECUTOR=${crosscompiling_emulator}"
                        -D "TEST_WORKING_DIR=${_WORKING_DIRECTORY}"
                        -D "TEST_EXTRA_ARGS=${_EXTRA_ARGS}"
                        -D "TEST_PROPERTIES=${_PROPERTIES}"
                        -D "TEST_PREFIX=${_TEST_PREFIX}"
                        -D "TEST_SUFFIX=${_TEST_SUFFIX}"
                        -D "TEST_FILTER=${_TEST_FILTER}"
                        -D "TEST_LIST={_TEST_LIST}"
                        -D "CTEST_FILE=${ctest_tests_file}"
                        -D "TEST_DISCOVERY_TIMEOUT=${_DISCOVERY_TIMEOUT}"
                        -D "TEST_XML_OUTPUT_DIR=${_XML_OUTPUT_DIR}"
                        -P "${_UNITTESTPP_DISCOVER_TESTS_SCRIPT}"
                VERBATIM
        )
        file(WRITE "${ctest_include_file}"
                "if(EXISTS \"${ctest_tests_file}\")\n"
                "    include(\"${ctest_tests_file}\")\n"
                "else()\n"
                "    add_test(${TARGET}_NOT_BUILT ${TARGET}_NOT_BUILT)\n"
                "endif()\n"
                )
    else()
        message(FATAL_ERROR "only POST_BUILD discovery mode supported")
    endif()

    set_property(DIRECTORY
            APPEND PROPERTY TEST_INCLUDE_FILES "${ctest_include_file}")

endfunction()


set(_UNITTESTPP_DISCOVER_TESTS_SCRIPT
        ${CMAKE_CURRENT_LIST_DIR}/UnitTestPP_DiscoverTests.cmake)


cmake_policy(POP)
