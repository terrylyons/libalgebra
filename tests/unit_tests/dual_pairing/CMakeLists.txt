


add_library(Libalgebra_tests_dual_pairing_tests OBJECT EXCLUDE_FROM_ALL
        test_tensor_dual_pairing.cpp
        )

target_link_libraries(Libalgebra_tests_dual_pairing_tests PUBLIC UnitTest++ Libalgebra::Libalgebra la_unittests::utilities)

add_library(la_unittests::dual_pairing_tests ALIAS Libalgebra_tests_dual_pairing_tests)

add_executable(dual_pairing_tests_exe)
target_link_libraries(dual_pairing_tests_exe PRIVATE la_unittests::main Libalgebra_tests_dual_pairing_tests)

add_test(NAME test_algebra
        COMMAND dual_pairing_tests_exe)