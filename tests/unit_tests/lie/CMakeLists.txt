

add_library(lie_tests OBJECT
        test_basis.cpp
        test_dense_lie_multiplication.cpp
        test_hall_basis.cpp
        test_lie_multiplication.cpp
        )

target_link_libraries(lie_tests PUBLIC UnitTest++ Libalgebra::Libalgebra la_unittests::utilities)

add_executable(lie_tests_exe "")
target_link_libraries(lie_tests_exe PRIVATE la_unittests::main lie_tests)
