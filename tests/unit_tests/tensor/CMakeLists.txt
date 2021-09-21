

add_library(tensor_tests OBJECT
        test_basis.cpp
        test_dense_tensor_functions.cpp
        test_dense_tensor_multiplication.cpp
        test_hybrid_tensor_functions.cpp
        test_hybrid_tensor_multiplication.cpp
        test_sparse_tensor_functions.cpp
        test_sparse_tensor_multiplication.cpp
        test_signature_calc.cpp
        test_tensor_functions.cpp
        test_tensor_multiplication.cpp
        test_tensor_size_info.cpp
)

target_link_libraries(tensor_tests PUBLIC UnitTest++ Libalgebra::Libalgebra la_unittests::utilities)



add_executable(tensor_tests_exe "")
target_link_libraries(tensor_tests_exe PUBLIC la_unittests::main tensor_tests)