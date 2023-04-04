add_executable(multiplication_bm
        multiplication_bm.cpp)

target_link_libraries(multiplication_bm PUBLIC Libalgebra::Libalgebra)

target_compile_options(multiplication_bm PRIVATE
        -g
        -pg
        -fno-omit-frame-pointer
        -march=native
#        -fopt-info-vec-missed
#        -ftree-vectorizer-verbose=2
#        -qopt-report-stdout
#        -qopt-report=max
        )

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "IntelLLVM")
#    target_compile_options(multiplication_bm PRIVATE
#            -qopt-report-stdout
#            -qopt-report=max)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
#    target_compile_options(multiplication_bm PRIVATE
#            -fopt-info
#            -fno-tree-loop-distribute-patterns
#            -ftree-vectorizer-verbose=2)
endif()


target_link_libraries(multiplication_bm PRIVATE
        benchmark::benchmark
        benchmark::benchmark_main
        )


find_package(Boost REQUIRED COMPONENTS program_options timer)



function(add_single_tm_benchmark width depth tile_letters cache_letters)
    set(fun_name multiplication_${width}_${depth}_${tile_letters}_${cache_letters}_single)
    add_executable(${fun_name} multiplication_single.cpp)


    target_link_libraries(${fun_name} PRIVATE
            benchmark::benchmark
            Libalgebra::Libalgebra
            Boost::program_options
            Boost::timer
            Boost::headers
            )
    target_compile_definitions(${fun_name} PRIVATE
            LA_TM_BM_WIDTH=${width}
            LA_TM_BM_DEPTH=${depth}
            LA_TM_BM_TILE_LETTERS=${tile_letters}
            LA_TM_BM_CACHE_LETTERS=${cache_letters})
    target_compile_options(${fun_name} PRIVATE
            -g -pg -fno-omit-frame-pointer)
endfunction()



add_single_tm_benchmark(4 13 3 0)
add_single_tm_benchmark(4 13 3 1)
add_single_tm_benchmark(4 13 3 2)
add_single_tm_benchmark(4 13 3 3)
