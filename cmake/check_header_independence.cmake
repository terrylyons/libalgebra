

macro(check_header_independence resultvar hpath target)

    string(REPLACE ".h" ".cpp" hpath_cpp ${hpath})

    file(WRITE "${CMAKE_BINARY_DIR}/test_compile/${hpath_cpp}" "
    #include <${hpath}>

    int main() {
        return 0;
    }
    ")

    if (WIN32)
        string(REPLACE "\\" "_" hpath_varname ${hpath})
    else ()
        string(REPLACE "/" "_" hpath_varname ${hpath})
    endif ()
    string(REPLACE ".h" "" hpath_varname ${hpath_varname})

    set(outvar test_compile_${hpath_varname})

    add_executable(${outvar} ${CMAKE_BINARY_DIR}/test_compile/${hpath_cpp})
    target_link_libraries(${outvar} PRIVATE ${target})
    target_include_directories(${outvar} PRIVATE ${CMAKE_SOURCE_DIR})

    #try_compile(${hpath_varname}_success
    #        ${CMAKE_BINARY_DIR}
    #        ${CMAKE_BINARY_DIR}/test_compile/${hpath_cpp}
    #        OUTPUT_VARIABLE ${hpath_varname}_output)
    #if (${hpath_varname}_success)
    #    message(STATUS "${hpath} - success")
    #else ()
    #    message(STATUS "${hpath} - failure")
    #    message(STATUS "${${hpath_varname}_output}")
    #endif ()

endmacro()


macro(check_all_header_independent resultvar path target)

    file(GLOB_RECURSE cahi_header_files RELATIVE ${path}/.. "${path}/*.h")
    foreach (cahi_p ${cahi_header_files})
        check_header_independence(resultvar ${cahi_p} target)
    endforeach ()

endmacro()