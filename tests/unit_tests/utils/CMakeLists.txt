

add_library(utils_tests OBJECT
        test_integer_maths.cpp
        )

target_link_libraries(utils_tests PUBLIC UnitTest++ Libalgebra::Libalgebra la_unittests::utilities)


add_executable(utils_tests_exe run_tests.cpp)
target_link_libraries(utils_tests_exe PRIVATE utils_tests)

