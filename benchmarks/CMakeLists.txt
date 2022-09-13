# See https://github.com/google/benchmark for instructions on how to install benchmark.
find_package(benchmark REQUIRED)

add_executable(antipode_bm
        antipode_bm.cpp
)

add_executable(multiplication_bm
        multiplication_bm.cpp
        )

target_link_libraries(antipode_bm PRIVATE
        benchmark::benchmark
        benchmark::benchmark_main
        )

target_link_libraries(multiplication_bm PRIVATE
        benchmark::benchmark
        benchmark::benchmark_main
        )