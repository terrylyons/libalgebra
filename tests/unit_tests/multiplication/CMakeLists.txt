

add_library(Libalgebra_tests_multiplication_tests OBJECT EXCLUDE_FROM_ALL
        test_multiplication.cpp
        )

target_link_libraries(Libalgebra_tests_multiplication_tests PUBLIC UnitTest++ Libalgebra::Libalgebra la_unittests::utilities)

add_library(la_unittests::multiplication_tests ALIAS Libalgebra_tests_multiplication_tests)

add_executable(multiplication_tests_exe)
target_link_libraries(multiplication_tests_exe PRIVATE la_unittests::main Libalgebra_tests_multiplication_tests)

add_test(NAME multiplication_tests COMMAND multiplication_tests_exe)