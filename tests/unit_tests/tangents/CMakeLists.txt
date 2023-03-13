


#add_library(Libalgebra_tests_tangents_tests OBJECT EXCLUDE_FROM_ALL
#        test_tangents.h
#        test_tangents.cpp
#        )
#
#target_link_libraries(Libalgebra_tests_tangents_tests
#        PUBLIC
#            UnitTest++
#            Libalgebra::Libalgebra
#            la_unittests::utilities
#        )
#
#add_library(la_unittests::tangents_tests ALIAS Libalgebra_tests_tangents_tests)
#
#
#add_executable(tangents_tests_exe)
#target_link_libraries(tangents_tests_exe
#        PUBLIC
#            la_unittests::main
#            Libalgebra_tests_tangents_tests
#        )
#
#
#add_test(NAME tangents_tests COMMAND tangents_tests_exe)
libalgebra_test(NAME tangents
        SOURCES
        test_tangents.cpp
        test_tangents.h
        )
