




add_library(Libalgebra_tests_hall_set_tests OBJECT EXCLUDE_FROM_ALL
        test_hall_set.cpp
        )

target_link_libraries(Libalgebra_tests_hall_set_tests PUBLIC UnitTest++ Libalgebra::Libalgebra la_unittests::utilities)

add_library(la_unittests::hall_set_tests ALIAS Libalgebra_tests_hall_set_tests)

add_executable(hall_set_tests_exe)
target_link_libraries(hall_set_tests_exe PRIVATE la_unittests::main Libalgebra_tests_hall_set_tests)

add_test(NAME hall_set_tests
        COMMAND hall_set_tests_exe)