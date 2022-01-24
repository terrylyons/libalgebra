

add_library(
        Libalgebra_tests_coefficient_tests OBJECT EXCLUDE_FROM_ALL
)

target_sources(Libalgebra_tests_coefficient_tests PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/test_float_field.cpp

        )

target_link_libraries(Libalgebra_tests_coefficient_tests PUBLIC
        UnitTest++
        Libalgebra::Libalgebra
        la_unittests::utilities)

add_library(la_unittests::coefficient_tests ALIAS Libalgebra_tests_coefficient_tests)

add_executable(coefficient_tests_exe)
target_link_libraries(coefficient_tests_exe PRIVATE
        la_unittests::main
        Libalgebra_tests_coefficient_tests)

add_test(NAME test_coefficient
        COMMAND coefficient_tests_exe)