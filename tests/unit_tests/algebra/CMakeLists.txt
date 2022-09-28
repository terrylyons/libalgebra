
libalgebra_test(NAME algebra
        SOURCES
        test_algebra_product.cpp
        )

#add_library(Libalgebra_tests_algebra_tests OBJECT EXCLUDE_FROM_ALL
#        test_algebra_product.cpp
#        )
#
#target_link_libraries(Libalgebra_tests_algebra_tests PUBLIC UnitTest++ Libalgebra::Libalgebra la_unittests::utilities)
#
#add_library(la_unittests::algebra_tests ALIAS Libalgebra_tests_algebra_tests)
#
#add_executable(algebra_tests_exe)
#target_link_libraries(algebra_tests_exe PRIVATE la_unittests::main Libalgebra_tests_algebra_tests)
#
#add_test(NAME test_algebra
#        COMMAND algebra_tests_exe)
