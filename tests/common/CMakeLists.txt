

add_library(libalgebra_tests_utilities STATIC
        reporter.cpp
        memfile.cpp
        reporter.h
        compat.h
        rng.h
        simple_basis.h
        time_and_details.h
        memfile.h
        random_vector_generator.h
        SHOW.h
        random_coeffs.h
        multi_test.h)

set_target_properties(libalgebra_tests_utilities PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(libalgebra_tests_utilities PUBLIC "${CMAKE_CURRENT_LIST_DIR}")

target_link_libraries(libalgebra_tests_utilities PUBLIC
        Libalgebra::Libalgebra
        UnitTest++::UnitTest++
        Boost::boost
        Boost::filesystem
        Boost::iostreams
        Boost::random
        OpenMP::OpenMP_CXX
        )

# Set compiler flags for reporting warnings when compiling tests
target_compile_options(libalgebra_tests_utilities PUBLIC
        $<$<CXX_COMPILER_ID:MSVC>: /W3>
        $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:
        -Wall
        -Wextra
        -Wpedantic
        -Wno-unused-parameter
        -Wno-unknown-pragmas
        >
        )


add_library(la_unittests::utilities ALIAS libalgebra_tests_utilities)


add_library(main STATIC run_tests.cpp)
target_link_libraries(main PUBLIC
        UnitTest++::UnitTest++
        libalgebra_tests_utilities)

add_library(la_unittests::main ALIAS main)
