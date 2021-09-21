
find_package(UnitTest++ CONFIG REQUIRED)
find_path(UnitTest++_INCLUDE_DIRS "UnitTest++/UnitTest++.h")
# include_directories(SYSTEM "${UnitTest++_INCLUDE_DIRS})
# target_link_libraries(main PRIVATE UnitTest++)

message(STATUS ">>> UnitTest++_INCLUDE_DIRS = ${UnitTest++_INCLUDE_DIRS}")

add_library(reporter STATIC
        ../common/reporter.cpp ../common/reporter.h "${UnitTest++_INCLUDE_DIRS}/UnitTest++/UnitTest++.h")
target_include_directories(reporter SYSTEM PRIVATE "${UnitTest++_INCLUDE_DIRS}")
target_link_libraries(reporter PUBLIC UnitTest++)


add_subdirectory(../unit_tests/algebra)
add_subdirectory(../unit_tests/lie)
add_subdirectory(../unit_tests/signature_tests)
add_subdirectory(../unit_tests/tensor)
add_subdirectory(../unit_tests/utils)
add_subdirectory(../unit_tests/vector)


add_library(integration_tests OBJECT
        ../integration_tests/AlgebraFunctionsTests.cpp
        ../integration_tests/HallSetTests.cpp
        ../integration_tests/LatticePathTests.cpp
        ../integration_tests/makebm.cpp
        ../common/memfile.cpp
        ../integration_tests/OMPSigsTests.cpp
        ../integration_tests/speed_tests.cpp
        ../integration_tests/tests_libalgebra-demo.cpp
        ../integration_tests/tests_TENSOR_LIE_CBH_MAPS.cpp
        ../integration_tests/TreeBufferHelper.cpp
        ../integration_tests/x64sigs.cpp
        )

target_link_libraries(integration_tests PUBLIC UnitTest++ reporter Libalgebra::Libalgebra)
target_link_libraries(integration_tests PRIVATE
        Boost::boost
        Boost::filesystem
        Boost::iostreams
        Boost::random
        OpenMP::OpenMP_CXX
        )

add_executable(integration_tests_exe ../integration_tests/LibAlgebraUnitTests.cpp)
target_link_libraries(integration_tests_exe PRIVATE integration_tests)