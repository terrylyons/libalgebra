

add_library(shuffle_tensor_tests OBJECT EXCLUDE_FROM_ALL
        test_shuffle_tensor.cpp
        )

target_link_libraries(shuffle_tensor_tests PUBLIC UnitTest++ Libalgebra::Libalgebra la_unittests::utilities)

add_executable(shuffle_tensor_tests_exe EXCLUDE_FROM_ALL)
target_link_libraries(shuffle_tensor_tests_exe PRIVATE la_unittests::main shuffle_tensor_tests)

add_test(NAME shuffle_tensor_tests COMMAND shuffle_tensor_tests_exe)