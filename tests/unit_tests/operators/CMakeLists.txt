

add_library(Libalgebra_tests_operators_tests OBJECT EXCLUDE_FROM_ALL

        )

target_link_libraries(Libalgebra_tests_operators_tests PUBLIC UnitTest++ Libalgebra::Libalgebra la_unittests::utilities)

add_library(la_unittests::operators_tests ALIAS Libalgebra_tests_operators_tests)

add_executable(operators_tests_exe)
target_link_libraries(operators_tests_exe PUBLIC la_unittests::main Libalgebra_tests_operators_tests)


add_test(NAME operators_tests COMMAND operators_tests_exe)