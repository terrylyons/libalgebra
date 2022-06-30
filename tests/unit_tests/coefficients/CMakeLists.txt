

add_library(
        Libalgebra_tests_coefficient_tests OBJECT EXCLUDE_FROM_ALL
)

target_sources(Libalgebra_tests_coefficient_tests PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/test_float_field.cpp
        ${CMAKE_CURRENT_LIST_DIR}/standard_coefficient_test_suite.h
        ${CMAKE_CURRENT_LIST_DIR}/unital_coefficient_test_suite.h
        ${CMAKE_CURRENT_LIST_DIR}/test_double_field.cpp
        ${CMAKE_CURRENT_LIST_DIR}/test_mpfloat_field.cpp
        ${CMAKE_CURRENT_LIST_DIR}/test_mpfloat_fallback_field.cpp
        # Annoyingly, there is a problem randomly generating rationals at the moment
        # ${CMAKE_CURRENT_LIST_DIR}/test_rational_field.cpp
        #${CMAKE_CURRENT_LIST_DIR}/test_rational_fallback_field.cpp
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

add_test(NAME coefficient_tests
        COMMAND coefficient_tests_exe)