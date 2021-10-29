

add_library(Libalgebra_tests_lie_tests OBJECT EXCLUDE_FROM_ALL
        test_basis.cpp
        test_dense_lie_multiplication.cpp
        test_hall_basis.cpp
        test_lie_multiplication.cpp
        "test_function_extension.cpp")

target_link_libraries(Libalgebra_tests_lie_tests PUBLIC UnitTest++ Libalgebra::Libalgebra la_unittests::utilities)

add_library(la_unittests::lie_tests ALIAS Libalgebra_tests_lie_tests)

add_executable(lie_tests_exe)
target_link_libraries(lie_tests_exe PRIVATE la_unittests::main Libalgebra_tests_lie_tests)

add_test(NAME lie_tests COMMAND lie_tests_exe)