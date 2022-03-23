

add_library(Libalgebra_tests_involute_tests OBJECT EXCLUDE_FROM_ALL
        test_involute.cpp
        )

target_link_libraries(Libalgebra_tests_involute_tests PUBLIC UnitTest++ Libalgebra::Libalgebra la_unittests::utilities)

add_library(la_unittests::involute_tests ALIAS Libalgebra_tests_involute_tests)

add_executable(involute_tests_exe)
target_link_libraries(involute_tests_exe PRIVATE la_unittests::main Libalgebra_tests_involute_tests)

add_test(NAME involute_tests COMMAND involute_tests_exe)