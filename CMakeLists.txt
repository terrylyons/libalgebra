cmake_minimum_required(VERSION 3.0)

# I took this from CATCH2's (V2.0) CMakeList.txt because
# this seems to be close to our use case here.
# https://github.com/catchorg/Catch2/blob/v2.x/CMakeLists.txt
if(NOT DEFINED PROJECT_NAME)
    set(NOT_SUBPROJECT ON)
else()
    set(NOT_SUBPROJECT OFF)
endif()


project(Libalgebra VERSION 1.0.0)


include(GNUInstallDirs)


set(CMAKE_CXX_STANDARD 98)

## we have a limited binary boost dependency
find_package(Boost REQUIRED COMPONENTS thread container)
## we also require the headers only boost-multiprecision so get path and check it is there
find_path(BOOST_MULTIPRECISION_INCLUDE_DIRS "boost/multiprecision/complex128.hpp")
# target_include_directories(main PRIVATE ${BOOST_MULTIPRECISION_INCLUDE_DIRS})


## we require an arbitrary precision arithmetic package
if(WIN32)
find_library(Bignum_LIBRARY mpir)
else()
find_library(Bignum_LIBRARY gmp)
endif()
find_path(Bignum_INCLUDE_DIR NAMES gmp.h gmpxx.h)
# include_directories(SYSTEM "${Bignum_INCLUDE_DIR}")

message(STATUS ">>> Bignum library: ${Bignum_LIBRARY}")
message(STATUS ">>> Bignum INCLUDE_DIRECTORIES: ${Bignum_INCLUDE_DIR}")

add_library(Libalgebra INTERFACE)

target_include_directories(Libalgebra
        INTERFACE
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        )

add_library(Libalgebra::Libalgebra ALIAS Libalgebra)

target_link_libraries(Libalgebra INTERFACE Boost::boost Boost::thread Boost::container "${Bignum_LIBRARY}")

if(NOT_SUBPROJECT)
    install(TARGETS Libalgebra
            EXPORT LibalgebraTargets
            INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
            )

    include(CMakePackageConfigHelpers)
    set(LIBALGEBRA_CMAKE_CONFIG_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Libalgebra")

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
    install(FILES   ${CMAKE_CURRENT_BINARY_DIR}/LibalgebraConfig.cmake
                    ${CMAKE_CURRENT_BINARY_DIR}/LibalgebraConfigVersion.cmake
            DESTINATION ${LIBALGEBRA_CMAKE_CONFIG_DESTINATION}
            )


    install(EXPORT LibalgebraTargets
            FILE LibalgebraTargets.cmake
            NAMESPACE Libalgebra::
            DESTINATION lib/cmake/Libalgebra
            )
endif()