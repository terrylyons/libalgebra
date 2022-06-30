

add_library(Libalgebra_tests_multipoly_tests OBJECT EXCLUDE_FROM_ALL)

target_link_libraries(Libalgebra_tests_multipoly_tests PUBLIC
        UnitTest++
        Libalgebra::Libalgebra
        la_unittests::utilities
        )
add_library(la_unittests::multipoly_tests ALIAS Libalgebra_tests_multipoly_tests)

target_sources(Libalgebra_tests_multipoly_tests PRIVATE
        test_multipoly_prod.cpp)

add_executable(multipoly_tests_exe)
target_link_libraries(multipoly_tests_exe PUBLIC la_unittests::main
        Libalgebra_tests_multipoly_tests)

add_test(NAME multipoly_tests COMMAND multipoly_tests_exe)
