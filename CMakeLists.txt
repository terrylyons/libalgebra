cmake_minimum_required(VERSION 3.1)

# I took this from CATCH2's (V2.0) CMakeList.txt because
# this seems to be close to our use case here.
# https://github.com/catchorg/Catch2/blob/v2.x/CMakeLists.txt
if(NOT DEFINED PROJECT_NAME)
    set(NOT_SUBPROJECT ON)
else()
    set(NOT_SUBPROJECT OFF)
endif()


project(Libalgebra VERSION 1.0.0)

option(LIBALGEBA_TESTING "Enable building test targets" OFF)


add_compile_options(
        "-Wpedantic"
        "-Wconversion"
        "-Wsign-compare"
)


include(GNUInstallDirs)
set(LIBALGEBRA_CMAKE_CONFIG_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Libalgebra")

#set(CMAKE_CXX_STANDARD 98)

## we have a limited binary boost dependency
find_package(Boost REQUIRED COMPONENTS thread container)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules/)
find_package(Bignum REQUIRED)



add_library(Libalgebra INTERFACE)

add_custom_target(Libalgebra.headers SOURCES
        libalgebra/basis/tags.h
        libalgebra/basis/basis.h
        libalgebra/coefficients/coefficients.h
        libalgebra/utils/integer_maths.h
        libalgebra/utils/iterators.h
        libalgebra/utils/meta.h
        libalgebra/utils/order_trait.h
        libalgebra/vectors/base_vector.h
        libalgebra/vectors/dense_storage.h
        libalgebra/vectors/dense_vector.h
        libalgebra/vectors/hybrid_vector.h
        libalgebra/vectors/iterators.h
        libalgebra/vectors/sparse_vector.h
        libalgebra/vectors/vector.h
        libalgebra/vectors/vectors.h
        libalgebra/_tensor_basis.h
        libalgebra/alg_types.h
        libalgebra/algebra.h
        libalgebra/basis_traits.h
        libalgebra/compat.h
        libalgebra/constlog2.h
        libalgebra/constpower.h
        libalgebra/hall_set.h
        libalgebra/implementation_types.h
        libalgebra/libalgebra.h
        libalgebra/lie.h
        libalgebra/lie_basis.h
        libalgebra/monomial_basis.h
        libalgebra/multi_polynomial.h
        libalgebra/multiplication_helpers.h
        libalgebra/poly_basis.h
        libalgebra/poly_lie.h
        libalgebra/poly_lie_basis.h
        libalgebra/polynomials.h
        libalgebra/tensor.h
        libalgebra/tensor_basis.h
        libalgebra/utils.h
        )

target_include_directories(Libalgebra
        INTERFACE
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        )
target_link_libraries(Libalgebra INTERFACE Boost::boost Boost::thread Boost::container Bignum::Bignum)



add_library(Libalgebra::Libalgebra ALIAS Libalgebra)


if(NOT_SUBPROJECT)
    install(TARGETS Libalgebra
            EXPORT LibalgebraTargets
            PUBLIC_HEADER
            INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
            )

    include(CMakePackageConfigHelpers)

    get_target_property(LIBALGEBRA_HEADER_SOURCES Libalgebra.headers SOURCES)
    install(DIRECTORY ${CMAKE_SOURCE_DIR}/libalgebra/
            DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/libalgebra)


    message(STATUS "Writing cmake config file")
    configure_package_config_file(
            LibalgebraConfig.cmake.in
            ${CMAKE_CURRENT_BINARY_DIR}/LibalgebraConfig.cmake
            INSTALL_DESTINATION ${LIBALGEBRA_CMAKE_CONFIG_DESTINATION}
            PATH_VARS CMAKE_INSTALL_INCLUDEDIR
    )

    message(STATUS "Writing cmake config version file")
    write_basic_package_version_file(
            LibalgebraConfigVersion.cmake
            VERSION ${PACKAGE_VERSION}
            COMPATIBILITY SameMajorVersion
            ARCH_INDEPENDENT
    )

    message(STATUS "Installing config files to ${LIBALGEBRA_CMAKE_CONFIG_DESTINATION}")


    install(EXPORT LibalgebraTargets
            FILE LibalgebraTargets.cmake
            NAMESPACE Libalgebra::
            DESTINATION ${LIBALGEBRA_CMAKE_CONFIG_DESTINATION}
            )

    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/LibalgebraConfig.cmake
            ${CMAKE_CURRENT_BINARY_DIR}/LibalgebraConfigVersion.cmake
            ${CMAKE_SOURCE_DIR}/cmake/Modules/FindBignum.cmake
            DESTINATION ${LIBALGEBRA_CMAKE_CONFIG_DESTINATION}
            )


endif()

if (LIBALGEBRA_TESTING)
    enable_testing()
    add_subdirectory(tests)
endif()