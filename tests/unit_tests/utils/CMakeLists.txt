

add_library(Libalgebra_tests_utils_tests OBJECT EXCLUDE_FROM_ALL
        test_integer_maths.cpp
        )

target_link_libraries(Libalgebra_tests_utils_tests PUBLIC UnitTest++ Libalgebra::Libalgebra la_unittests::utilities)

add_library(la_unittests::utils_tests ALIAS Libalgebra_tests_utils_tests)

add_executable(utils_tests_exe)
target_link_libraries(utils_tests_exe PRIVATE la_unittests::main Libalgebra_tests_utils_tests)

add_test(NAME utils_tests COMMAND utils_tests_exe)