

add_library(utilities STATIC EXCLUDE_FROM_ALL
        reporter.cpp
        memfile.cpp
        reporter.h
        compat.h
        rng.h
        simple_basis.h
        small_vector.h
        time_and_details.h
        memfile.h
        )

target_include_directories(utilities PUBLIC "${CMAKE_CURRENT_LIST_DIR}")

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


add_library(main STATIC EXCLUDE_FROM_ALL run_tests.cpp)
target_link_libraries(main PUBLIC utilities)

add_library(la_unittests::main ALIAS main)