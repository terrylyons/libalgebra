add_executable(antipode_bm
        antipode_bm.cpp
        )

target_link_libraries(antipode_bm PUBLIC Libalgebra::Libalgebra)

#if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
target_compile_options(antipode_bm PRIVATE
        -fno-omit-frame-pointer
        -g -pg
        -march=native)
#endif()

target_link_libraries(antipode_bm PRIVATE
        benchmark::benchmark
        benchmark::benchmark_main
        )
