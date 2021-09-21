

add_library(utilities STATIC
        reporter.cpp
        memfile.cpp
        )

target_sources(utilities PUBLIC
        reporter.h
        compat.h
        rng.h
        simple_basis.h
        small_vector.h
        time_and_details.h
        memfile.h
        )

target_link_libraries(utilities PUBLIC
        Libalgebra::Libalgebra
        UnitTest++
        Boost::boost
        Boost::filesystem
        Boost::iostreams
        Boost::random
        OpenMP::OpenMP_CXX
        )

add_library(la_unittests::utilities ALIAS utilities)