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

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# This helps IDEs that can arrange projects into folders (apparently)
# See https://cmake.org/cmake/help/v3.22/prop_gbl/USE_FOLDERS.html
# from https://github.com/Lectem/cpp-boilerplate/blob/master/CMakeLists.txt
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

include(vcpkg OPTIONAL)
include(GNUInstallDirs)

set(LIBALGEBRA_CMAKE_CONFIG_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Libalgebra")


## we have a limited binary boost dependency
set(Boost_NO_WARN_NEW_VERSIONS ON)
find_package(Boost REQUIRED COMPONENTS thread)

if (NOT LIBALGEBRA_NO_SERIALIZATION)
    find_package(Boost OPTIONAL_COMPONENTS serialization)
endif()


list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/)
find_package(Bignum)

# Interface libraries sometimes complain if you list sources, see below.
add_library(Libalgebra INTERFACE)

# For IDEs that use target sources lists to generate project information
# we add a custom target that attaches these sources.
add_custom_target(Libalgebra_headers SOURCES
        libalgebra/basis/tags.h
        libalgebra/basis/basis.h
        libalgebra/basis/key_iterators.h
        libalgebra/coefficients/coefficients.h
        libalgebra/utils/caching_tags.h
        libalgebra/utils/function_extension_cache_base.h
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
        libalgebra/operators/operators.h
        libalgebra/operators/dot_product_implementations.h
        libalgebra/operators/functionals.h
        libalgebra/operators/free_extension.h
        libalgebra/operators/tensor_operator.h
        libalgebra/operators/lie_inner_product.h
        libalgebra/operators/scalar_multiply_operator.h
        libalgebra/operators/sum_operator.h
        libalgebra/operators/composition_operator.h
        libalgebra/operators/multi_linear_operators.h
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
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
        )
target_link_libraries(Libalgebra INTERFACE
        Boost::boost
        Boost::thread
        )

if (TARGET Boost::serialization)
    target_link_libraries(Libalgebra INTERFACE Boost::serialization)
    target_compile_definitions(Libalgebra INTERFACE LIBALGEBRA_ENABLE_SERIALIZATION)
endif()

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

if (LIBALGEBRA_TESTING)
    # Set the C++ standard to c++11 required, so we get better linting in IDEs
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)

    enable_testing()
    add_subdirectory(tests)
endif ()