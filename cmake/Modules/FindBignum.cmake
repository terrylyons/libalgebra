
# Basic find module for bignum library.
# I'm using https://gitlab.kitware.com/cmake/community/-/wikis/doc/tutorials/How-To-Find-Libraries
# as reference here.


if (WIN32)
    find_library(Bignum_LIBRARY mpir)
    find_path(Bignum_INCLUDE_DIR gmp.h gmpxx.h)

else()
    find_package(PkgConfig)
    pkg_check_modules(PC_GMP QUIET gmp)

    set(GMP_DEFINITIONS ${PC_GMP_CFLAGS_OTHER})

    find_library(Bignum_LIBRARY NAMES gmp
            HINTS ${PC_GMP_LIBDIR} ${PC_GMP_LIBRARY_DIRS})

    find_path(Bignum_INCLUDE_DIR NAMES gmp.h gmpxx.h
            HINTS ${PC_GMP_INCLUDEDIR} ${PC_GMP_INCLUDE_DIRS}
            )
endif()


include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Bignum DEFAULT_MSG Bignum_LIBRARY Bignum_INCLUDE_DIR)
mark_as_advanced(Bignum_LIBRARY Bignum_INCLUDE_DIR)

set(Bignum_LIBRARIES ${Bignum_LIBRARY})
set(Bignum_INCLUDE_DIRS ${Bignum_INCLUDE_DIR})

if (Bignum_FOUND AND NOT TARGET Bignum::Bignum)
    # Add an imported target that carries all the information we need. We can just link this target in other projects
    add_library(Bignum::Bignum SHARED IMPORTED GLOBAL)
    set_property(TARGET Bignum::Bignum PROPERTY IMPORTED_LOCATION "${Bignum_LIBRARY}")
    target_include_directories(Bignum::Bignum INTERFACE ${Bignum_INCLUDE_DIR})
    target_compile_definitions(Bignum::Bignum INTERFACE ${GMP_DEFINITIONS})

    # MSVC requires this for some reason?
    set_property(TARGET Bignum::Bignum PROPERTY IMPORTED_IMPLIB "${Bignum_LIBRARY}")
endif()