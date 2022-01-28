

add_library(Libalgebra_tests_integration_tests OBJECT EXCLUDE_FROM_ALL
        alg_framework.h
        brown_path_increments.h
        log2ceil.h
        makebm.h
        SigHelpers.h
        TreeBufferHelper.h
        AlgebraFunctionsTests.cpp
        HallSetTests.cpp
        LatticePathTests.cpp
        makebm.cpp
        OMPSigsTests.cpp
        speed_tests.cpp
        tests_libalgebra-demo.cpp
        tests_TENSOR_LIE_CBH_MAPS.cpp
        TreeBufferHelper.cpp
        x64sigs.cpp
        )


target_link_libraries(Libalgebra_tests_integration_tests PUBLIC
        la_unittests::utilities
        Libalgebra::Libalgebra
        )
add_library(la_unittests::integration_tests ALIAS Libalgebra_tests_integration_tests)

add_executable(integration_tests_exe)
target_link_libraries(integration_tests_exe PRIVATE la_unittests::main Libalgebra_tests_integration_tests)

add_test(NAME integration_tests COMMAND integration_tests_exe)
add_subdirectory(signature_tests)