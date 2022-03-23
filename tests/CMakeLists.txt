cmake_minimum_required(VERSION 3.12)


if (WIN32)
    ## better remove this and add it to the command line as machine specific
    #list(APPEND CMAKE_MODULE_PATH  "c:/dev/vcpkg/scripts/buildsystems/")
    add_compile_options(/bigobj /fp:strict)
else ()
endif ()


find_package(UnitTest++ CONFIG REQUIRED)

# include_directories(SYSTEM "${UnitTest++_INCLUDE_DIRS})
# target_link_libraries(main PRIVATE UnitTest++)

find_package(Boost COMPONENTS thread filesystem iostreams chrono random container REQUIRED)

find_package(OpenMP REQUIRED)

option(LA_TESTING_INCLUDE_INTEGRATION OFF)


message(STATUS ">>> Boost_INCLUDE_DIRS = ${Boost_INCLUDE_DIRS}")
message(STATUS ">>> UnitTest++_INCLUDE_DIRS = ${UnitTest++_INCLUDE_DIRS}")
message(STATUS ">>> OpenMP_CXX_INCLUDE_DIRS = ${OpenMP_CXX_INCLUDE_DIRS}")
message(STATUS ">>> OpenMP_CXX_FLAGS = ${OpenMP_CXX_FLAGS}")
message(STATUS ">>> Boost libraries: ${Boost_LIBRARIES}")

add_subdirectory(common)
if (LA_TESTIN_INCLUDE_INTEGRATION)
    add_subdirectory(integration_tests)
endif ()

add_subdirectory(unit_tests)
add_executable(all_tests EXCLUDE_FROM_ALL)
target_link_libraries(all_tests PRIVATE
        la_unittests::main
        la_unittests::algebra_tests
        la_unittests::lie_tests
        la_unittests::tensor_tests
        la_unittests::utils_tests
        la_unittests::vector_tests
        la_unittests::hall_set_tests
        la_unittests::serialization_tests
        la_unittests::shuffle_tensor_tests
        la_unittests::operators_tests
        )

if (LA_TESTING_INCLUDE_INTEGRATION)

    target_link_libraries(all_tests PRIVATE
            la_unittests::integration_tests
            la_unittests::signature_tests
            )

endif ()