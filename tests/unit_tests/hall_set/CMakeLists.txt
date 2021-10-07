




add_library(hall_set_tests OBJECT EXCLUDE_FROM_ALL
        test_hall_set.cpp
        )

target_link_libraries(hall_set_tests PUBLIC UnitTest++ Libalgebra::Libalgebra la_unittests::utilities)

add_executable(hall_set_tests_exe EXCLUDE_FROM_ALL)
target_link_libraries(hall_set_tests_exe PRIVATE la_unittests::main hall_set_tests)

add_test(NAME test_algebra
        COMMAND hall_set_tests_exe)