
find_package(Boost COMPONENTS REQUIRED serialization)

libalgebra_test(NAME serialization
        SOURCES
        dense_vec_serialization.cpp
        fixture.h
        fixture.cpp
        temporary_directory.cpp
        temporary_directory.h
        test_dense_storage.cpp
        sparse_vector_serialization.cpp
        test_higher_classes.cpp
        test_hall_set_serialize.cpp
        hybrid_vector_serialization.cpp
        DEPS
        Boost::boost
        Boost::serialization
        DEFN
            LIBALGEBRA_ENABLE_SERIALIZATION
        )
