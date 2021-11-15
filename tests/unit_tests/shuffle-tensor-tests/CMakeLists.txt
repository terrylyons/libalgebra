

add_library(Libalgebra_tests_shuffle_tensor_tests OBJECT EXCLUDE_FROM_ALL
        test_shuffle_tensor.cpp
        )

target_link_libraries(Libalgebra_tests_shuffle_tensor_tests PUBLIC UnitTest++ Libalgebra::Libalgebra la_unittests::utilities)

add_library(la_unittests::shuffle_tensor_tests ALIAS Libalgebra_tests_shuffle_tensor_tests)

add_executable(shuffle_tensor_tests_exe)
target_link_libraries(shuffle_tensor_tests_exe PRIVATE la_unittests::main Libalgebra_tests_shuffle_tensor_tests)

add_test(NAME shuffle_tensor_tests COMMAND shuffle_tensor_tests_exe)