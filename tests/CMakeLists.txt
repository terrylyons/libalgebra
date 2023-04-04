cmake_minimum_required(VERSION 3.12)



if (WIN32)
    ## better remove this and add it to the command line as machine specific
    #list(APPEND CMAKE_MODULE_PATH  "c:/dev/vcpkg/scripts/buildsystems/")
    add_compile_options(/bigobj /fp:strict)
else ()
endif ()


find_package(UnitTest++ CONFIG REQUIRED)

if (NOT TARGET UnitTest++::UnitTest++ AND TARGET UnitTest++)
    add_library(UnitTest++::UnitTest++ ALIAS UnitTest++)
endif()

include(libalgebra_test_helper)
# include_directories(SYSTEM "${UnitTest++_INCLUDE_DIRS})
# target_link_libraries(main PRIVATE UnitTest++)

find_package(Boost REQUIRED COMPONENTS thread filesystem iostreams chrono random container )

find_package(OpenMP REQUIRED)

option(LA_TESTING_INCLUDE_INTEGRATION OFF)


message(STATUS ">>> Boost_INCLUDE_DIRS = ${Boost_INCLUDE_DIRS}")
message(STATUS ">>> UnitTest++_INCLUDE_DIRS = ${UnitTest++_INCLUDE_DIRS}")
message(STATUS ">>> OpenMP_CXX_INCLUDE_DIRS = ${OpenMP_CXX_INCLUDE_DIRS}")
message(STATUS ">>> OpenMP_CXX_FLAGS = ${OpenMP_CXX_FLAGS}")
message(STATUS ">>> Boost libraries: ${Boost_LIBRARIES}")

add_subdirectory(common)
if (LA_TESTING_INCLUDE_INTEGRATION)
    add_subdirectory(integration_tests)
endif ()

add_subdirectory(unit_tests)
#add_executable(all_tests common/run_tests.cpp)
