

add_library(vector_tests OBJECT
        test_dense.cpp
        test_dense_storage.cpp
        test_hybrid.cpp
        test_sparse.cpp
        )
target_link_libraries(vector_tests PUBLIC Libalgebra::Libalgebra UnitTest++ la_unittests::utilities)

add_executable(vector_tests_exe "")
target_link_libraries(vector_tests_exe PRIVATE la_unittests::main vector_tests)