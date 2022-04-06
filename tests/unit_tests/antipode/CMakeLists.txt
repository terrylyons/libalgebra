

add_library(Libalgebra_tests_antipode_tests OBJECT EXCLUDE_FROM_ALL
        test_antipode.cpp
        )

target_link_libraries(Libalgebra_tests_antipode_tests PUBLIC UnitTest++ Libalgebra::Libalgebra la_unittests::utilities)

add_library(la_unittests::antipode_tests ALIAS Libalgebra_tests_antipode_tests)

add_executable(antipode_tests_exe)
target_link_libraries(antipode_tests_exe PRIVATE la_unittests::main Libalgebra_tests_antipode_tests)

add_test(NAME antipode_tests COMMAND antipode_tests_exe)