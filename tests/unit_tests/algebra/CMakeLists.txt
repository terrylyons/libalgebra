


add_library(algebra_tests OBJECT
        test_algebra_product.cpp
        )

target_link_libraries(algebra_tests PUBLIC UnitTest++ Libalgebra::Libalgebra la_unittests::utilities)

add_executable(algebra_tests_exe run_tests.cpp)
target_link_libraries(algebra_tests_exe PRIVATE algebra_tests)

