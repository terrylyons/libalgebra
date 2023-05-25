cmake_minimum_required(VERSION 3.1)

# I took this from CATCH2's (V2.0) CMakeList.txt because
# this seems to be close to our use case here.
# https://github.com/catchorg/Catch2/blob/v2.x/CMakeLists.txt
if (NOT DEFINED PROJECT_NAME)
    set(NOT_SUBPROJECT ON)
else ()
    set(NOT_SUBPROJECT OFF)
endif ()


project(Libalgebra VERSION 2.0.0)

option(LIBALGEBRA_TESTING "Enable building test targets" OFF)
option(LIBALGEBRA_NO_SERIALIZATION "Turn off building serialization dependency" OFF)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# This helps IDEs that can arrange projects into folders (apparently)
# See https://cmake.org/cmake/help/v3.22/prop_gbl/USE_FOLDERS.html
# from https://github.com/Lectem/cpp-boilerplate/blob/master/CMakeLists.txt
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

include(vcpkg OPTIONAL)
include(GNUInstallDirs)

set(LIBALGEBRA_CMAKE_CONFIG_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Libalgebra")

# -DLIBALGEBRA_MAX_TILE_LETTERS=X
if (DEFINED LIBALGEBRA_MAX_TILE_LETTERS)
    message("You have set LIBALGEBRA_MAX_TILE_LETTERS to ${LIBALGEBRA_MAX_TILE_LETTERS}")
    add_definitions(-DLIBALGEBRA_MAX_TILE_LETTERS=${LIBALGEBRA_MAX_TILE_LETTERS})
endif ()

## we have a limited binary boost dependency
set(Boost_NO_WARN_NEW_VERSIONS ON)
find_package(Boost REQUIRED)

if (NOT LIBALGEBRA_NO_SERIALIZATION)
    find_package(Boost OPTIONAL_COMPONENTS serialization)
endif ()


list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)
find_package(Bignum)

find_package(OpenMP)

#message(STATUS ">>> OpenMP_CXX_INCLUDE_DIRS = ${OpenMP_CXX_INCLUDE_DIRS}")
#message(STATUS ">>> OpenMP_CXX_FLAGS = ${OpenMP_CXX_FLAGS}")

# Interface libraries sometimes complain if you list sources, see below.
add_library(Libalgebra INTERFACE)


# For IDEs that use target sources lists to generate project information
# we add a custom target that attaches these sources.
add_custom_target(Libalgebra_headers SOURCES
        libalgebra/tags.h
        libalgebra/basis.h
        libalgebra/key_iterators.h
        libalgebra/coefficients.h
        libalgebra/complex.h
        libalgebra/mpfloat_coefficients.h
        libalgebra/rational_coefficients.h
        libalgebra/detail/caching_tags.h
        libalgebra/detail/function_extension_cache_base.h
        libalgebra/detail/integer_maths.h
        libalgebra/detail/meta.h
        libalgebra/detail/order_trait.h
        libalgebra/base_vector.h
        libalgebra/dense_storage.h
        libalgebra/dense_vector.h
        libalgebra/hybrid_vector.h
        libalgebra/iterators.h
        libalgebra/sparse_vector.h
        libalgebra/vector.h
        libalgebra/vectors.h
        libalgebra/_tensor_basis.h
        libalgebra/alg_types.h
        libalgebra/algebra.h
        libalgebra/base_basis.h
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
        libalgebra/operators.h
        libalgebra/dot_product_implementations.h
        libalgebra/functionals.h
        libalgebra/free_extension.h
        libalgebra/tensor_operator.h
        libalgebra/lie_inner_product.h
        libalgebra/scalar_multiply_operator.h
        libalgebra/sum_operator.h
        libalgebra/composition_operator.h
        libalgebra/multi_linear_operators.h
        libalgebra/poly_basis.h
        libalgebra/poly_lie.h
        libalgebra/poly_lie_basis.h
        libalgebra/polynomials.h
        libalgebra/tensor.h
        libalgebra/reconfigure_operator.h
        libalgebra/tensor_basis.h
        libalgebra/utils.h
        libalgebra/alternative_multiplications.h
        libalgebra/area_tensor_basis.h
        libalgebra/area_tensor_multiplication.h
        libalgebra/half_shuffle_tensor_basis.h
        libalgebra/half_shuffle_tensor_multiplication.h
        libalgebra/vector_bundle.h
        libalgebra/detail/reversing_permutation.h
        )


target_include_directories(Libalgebra
        INTERFACE
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
        )
target_link_libraries(Libalgebra INTERFACE Boost::boost)

if (TARGET Boost::serialization)
    target_link_libraries(Libalgebra INTERFACE Boost::serialization)
    target_compile_definitions(Libalgebra INTERFACE LIBALGEBRA_ENABLE_SERIALIZATION)
endif ()

# Set the _DEBUG flag when we are building in Debug configuration.
# This macro should probably be renamed to avoid collisions with compiler
# automatic definitions
target_compile_definitions(Libalgebra INTERFACE $<$<CONFIG:Debug>:_DEBUG>)

set_property(TARGET Libalgebra PROPERTY CXX_STANDARD 11)
set_property(TARGET Libalgebra PROPERTY CXX_STANDARD_REQUIRED ON)
# Set the __cplusplus macro in MSVC to a current compiler value
if (MSVC)
    target_compile_options(Libalgebra INTERFACE "/Zc:__cplusplus")
endif ()
if (TARGET Bignum::Bignum)
    target_link_libraries(Libalgebra INTERFACE Bignum::Bignum)
else ()
    target_compile_definitions(Libalgebra INTERFACE LIBALGEBRA_NO_GMP)
endif ()


add_library(Libalgebra::Libalgebra ALIAS Libalgebra)

# Build the documentation using Doxygen
if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.9)
    find_package(Doxygen OPTIONAL_COMPONENTS dot mscgen dia)
    if (Doxygen_FOUND)
        message(STATUS "Doxygen found, building docs")
        set(DOXYGEN_USE_MDFILE_AS_MAINPAGE README.md)
        doxygen_add_docs(
                doc
                README.md libalgebra
                COMMENT "Generate man pages"
        )
    else ()
        message(STATUS "Doxygen not found")
    endif ()
endif ()


if (NOT_SUBPROJECT)
    install(TARGETS Libalgebra
            EXPORT LibalgebraTargets
            PUBLIC_HEADER
            INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
            )

    include(CMakePackageConfigHelpers)

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
endif ()

enable_testing()
if (LIBALGEBRA_TESTING)
    add_subdirectory(tests)
endif ()

if (LIBALGEBRA_BENCHMARKING)
    add_subdirectory(benchmarks)
endif ()
