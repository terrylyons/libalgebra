

add_library(integration_tests OBJECT
        alg_framework.h
        brown_path_increments.h
        log2ceil.h
        makebm.h
        SHOW.h
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


target_link_libraries(integration_tests PUBLIC
        la_unittests::utilities
        Libalgebra::Libalgebra
        )


add_executable(integration_tests_exe "")
target_link_libraries(integration_tests_exe PRIVATE la_unittests::main integration_tests)