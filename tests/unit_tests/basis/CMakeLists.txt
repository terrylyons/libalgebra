


add_library(Libalgebra_tests_basis_tests OBJECT EXCLUDE_FROM_ALL basis_iteration_suite.h test_tensor_basis_iteration.cpp test_lie_basis_iteration.cpp)

target_sources(Libalgebra_tests_basis_tests PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/test_basis_iteration.cpp
        )

target_link_libraries(Libalgebra_tests_basis_tests PUBLIC UnitTest++ Libalgebra::Libalgebra la_unittests::utilities)

add_library(la_unittests::basis_tests ALIAS Libalgebra_tests_basis_tests)


add_executable(basis_tests_exe)
target_link_libraries(basis_tests_exe PRIVATE la_unittests::main Libalgebra_tests_basis_tests)


add_test(NAME test_basis COMMAND basis_tests_exe)