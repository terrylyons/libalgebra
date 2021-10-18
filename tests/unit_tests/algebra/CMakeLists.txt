


add_library(algebra_tests OBJECT EXCLUDE_FROM_ALL
        test_algebra_product.cpp
        )

target_link_libraries(algebra_tests PUBLIC UnitTest++ Libalgebra::Libalgebra la_unittests::utilities)

add_executable(algebra_tests_exe)
target_link_libraries(algebra_tests_exe PRIVATE la_unittests::main algebra_tests)

add_test(NAME test_algebra
        COMMAND algebra_tests_exe)