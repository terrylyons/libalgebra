
include(${CMAKE_CURRENT_LIST_DIR}/Modules/UnitTestPP.cmake)

set(LIBALGEBRA_TEST_SUITE_LIBS "" CACHE INTERNAL "test libs")


function(libalgebra_test)

    cmake_parse_arguments(LA_TESTS
            ""                                               # options
            "NAME"                                           # one value args
            "SOURCES;DEPS;DEFN"  # multivalue args
            ${ARGN}
            )


    set(_lib_name "LA_TESTING_${LA_TESTS_NAME}")

    set(_srcs "")
    foreach(_file IN LISTS LA_TESTS_SOURCES)
        if (NOT IS_ABSOLUTE _file)
            set(_file "${CMAKE_CURRENT_LIST_DIR}/${_file}")
        endif()
        if (NOT EXISTS "${_file}")
            message(FATAL_ERROR "File ${_file} does not exist")
        endif()
        list(APPEND _srcs ${_file})
    endforeach()

    add_executable(${_lib_name} "${_srcs}")
#    add_library(${_lib_name} OBJECT "${_srcs}")
#    set(LIBALGEBRA_TEST_SUITE_LIBS ${LIBALGEBRA_TEST_SUITE_LIBS} "${_lib_name}" CACHE INTERNAL "test_libs")
    message(STATUS "Adding tests suite ${LA_TESTS_NAME}")

    set_target_properties(${_lib_name} PROPERTIES
            POSITION_INDEPENDENT_CODE ON
            LINKER_LANGUAGE CXX)

    target_link_libraries(${_lib_name} PUBLIC
            UnitTest++::UnitTest++
            Libalgebra::Libalgebra
            la_unittests::utilities
            la_unittests::main
            "${LA_TESTS_DEPS}"
            )
    if (LA_TESTS_DEFN)
        target_compile_definitions(${_lib_name} PUBLIC "${LA_TESTS_DEFN}")
    endif()

    UnitTestPP_discover_tests(${_lib_name})


endfunction()
