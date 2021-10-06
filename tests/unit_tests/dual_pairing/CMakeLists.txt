


add_library(dual_pairing_tests OBJECT EXCLUDE_FROM_ALL
        test_tensor_dual_pairing.cpp
        )

target_link_libraries(dual_pairing_tests PUBLIC UnitTest++ Libalgebra::Libalgebra la_unittests::utilities)

add_executable(dual_pairing_tests_exe EXCLUDE_FROM_ALL)
target_link_libraries(dual_pairing_tests_exe PRIVATE la_unittests::main dual_pairing_tests)

add_test(NAME test_algebra
        COMMAND dual_pairing_tests_exe)