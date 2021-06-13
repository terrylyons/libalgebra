cmake_minimum_required(VERSION 3.0)

project(Libalgebra VERSION 1.0.0)


include(GNUInstallDirs)


set(CMAKE_CXX_STANDARD 98)

find_package(Boost COMPONENTS thread)


add_library(Libalgebra INTERFACE ${CMAKE_CURRENT_LIST_DIR})

target_include_directories(Libalgebra
        INTERFACE
            $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/libalgebra>
        )

add_library(Libalgebra::Libalgebra ALIAS Libalgebra)

target_link_libraries(Libalgebra INTERFACE Boost)


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
