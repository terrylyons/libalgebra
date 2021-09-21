



add_library(signature_tests OBJECT
        test_brownian_path.cpp
        test_innovative_path.cpp
        )

target_link_libraries(signature_tests PUBLIC UnitTest++ Libalgebra::Libalgebra la_unittests::utilities)



add_executable(signature_tests_exe run_tests.cpp)
target_link_libraries(signature_tests_exe PRIVATE signature_tests)