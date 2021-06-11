cmake_minimum_required(VERSION 3.0)

project(libalgebra VERSION 1.0.0)





set(CMAKE_CXX_STANDARD 98)

find_package(Boost COMPONENTS thread)


add_library(Libalgebra INTERFACE
        libalgebra/libalgebra.h
        libalgebra/_tensor_basis.h
        libalgebra/algebra.h
        libalgebra/basis_traits.h
        libalgebra/constlog2.h
        libalgebra/constpower.h
        libalgebra/implimentation_types.h
        libalgebra/lie.h
        libalgebra/lie_basis.h
        libalgebra/monomial_basis.h
        libalgebra/multi_polynomial.h
        libalgebra/poly_basis.h
        libalgebra/poly_lie.h
        libalgebra/poly_lie_basis.h
        libalgebra/polynomials.h
        libalgebra/sparse_vector.h
        libalgebra/tensor.h
        libalgebra/tensor_basis.h
        libalgebra/utils.h
        )
target_link_libraries(Libalgebra INTERFACE Boost)

include(GNUInstallDirs)
install(TARGETS Libalgebra
        EXPORT LibalgebraTargets
        INCLUDES DESTINATION include
        )

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        LibalgebraConfigVersion.cmake
        VERSION ${PACKAGE_VERSION}
        COMPATIBILITY SameMajorVersion
)

install(EXPORT LibalgebraTargets
        FILE LibalgebraTargets.cmake
        NAMESPACE Libalgebra::
        DESTINATION lib/cmake/Libalgebra
        )

configure_file(LibalgebraConfig.cmake.in LibalgebraConfig.cmake @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/LibalgebraConfig.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/LibalgebraConfigVersion.cmake"
        DESTINATION lib/cmake/Libalgebra
        )