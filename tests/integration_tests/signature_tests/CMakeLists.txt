



add_library(Libalgebra_tests_signature_tests OBJECT EXCLUDE_FROM_ALL
        test_brownian_path.cpp
        test_innovative_path.cpp
        )

target_link_libraries(Libalgebra_tests_signature_tests PUBLIC UnitTest++ Libalgebra::Libalgebra la_unittests::utilities)

add_library(la_unittests::signature_tests ALIAS Libalgebra_tests_signature_tests)

add_executable(signature_tests_exe)
target_link_libraries(signature_tests_exe PRIVATE la_unittests::main Libalgebra_tests_signature_tests)

add_test(NAME signature_tests COMMAND signature_tests_exe)