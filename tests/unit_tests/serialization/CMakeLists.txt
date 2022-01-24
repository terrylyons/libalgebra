
find_package(Boost COMPONENTS REQUIRED serialization)

add_library(Libalgebra_tests_serialization_tests OBJECT EXCLUDE_FROM_ALL
        dense_vec_serialization.cpp
        fixture.h
        fixture.cpp
        temporary_directory.cpp
        temporary_directory.h
        test_dense_storage.cpp sparse_vector_serialization.cpp test_higher_classes.cpp test_hall_set_serialize.cpp hybrid_vector_serialization.cpp)

target_link_libraries(Libalgebra_tests_serialization_tests PRIVATE
        Boost::boost
        Boost::serialization
        Libalgebra::Libalgebra)

add_library(la_unittests::serialization_tests ALIAS Libalgebra_tests_serialization_tests)


add_executable(serialization_tests_exe)
target_link_libraries(serialization_tests_exe PRIVATE
        la_unittests::main
        la_unittests::utilities
        Libalgebra_tests_serialization_tests)


add_test(NAME serialization_tests COMMAND serialization_tests_exe)