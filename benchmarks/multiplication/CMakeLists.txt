add_executable(multiplication_bm
        multiplication_bm.cpp
        )

target_link_libraries(multiplication_bm PUBLIC Libalgebra::Libalgebra)

target_compile_options(multiplication_bm PRIVATE
        -g
        -pg
        -fno-omit-frame-pointer)

target_link_libraries(multiplication_bm PRIVATE
        benchmark::benchmark
        benchmark::benchmark_main
        )
