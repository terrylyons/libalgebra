
add_library(Libalgebra_tests_half_shuffle_tensor_tests OBJECT EXCLUDE_FROM_ALL
        test_half_shuffle.cpp
        )

target_link_libraries(Libalgebra_tests_half_shuffle_tensor_tests PUBLIC
        UnitTest++
        Libalgebra::Libalgebra
        la_unittests::utilities)

add_library(la_unittests::half_shuffle_tests ALIAS Libalgebra_tests_half_shuffle_tensor_tests)

add_executable(half_shuffle_tests_exe)

target_link_libraries(half_shuffle_tests_exe PRIVATE
        la_unittests::main
        Libalgebra_tests_half_shuffle_tensor_tests
        )


add_test(NAME half_shuffle_tests COMMAND half_shuffle_tests_exe)
