add_executable(antipode_bm
        antipode_bm.cpp
        )

target_link_libraries(antipode_bm PUBLIC Libalgebra::Libalgebra)

target_link_libraries(antipode_bm PRIVATE
        benchmark::benchmark
        benchmark::benchmark_main
        )