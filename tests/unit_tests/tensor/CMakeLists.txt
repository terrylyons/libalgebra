
libalgebra_test(NAME tensor
        SOURCES
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
        test_tensor_creation.cpp
        )
#
#
#add_library(Libalgebra_tests_tensor_tests OBJECT EXCLUDE_FROM_ALL

#        )
#
#target_link_libraries(Libalgebra_tests_tensor_tests PUBLIC UnitTest++ Libalgebra::Libalgebra la_unittests::utilities)
#
#add_library(la_unittests::tensor_tests ALIAS Libalgebra_tests_tensor_tests)
#
#add_executable(tensor_tests_exe)
#target_link_libraries(tensor_tests_exe PUBLIC la_unittests::main Libalgebra_tests_tensor_tests)
#
#
#add_test(NAME tensor_tests COMMAND tensor_tests_exe)
