
include(CheckCXXSourceCompiles)
include(FindPackageHandleStandardArgs)
include(FindPackageMessage)

set(SEQAN3_BENCHMARK_SRC_DIR "${PROJECT_BINARY_DIR}/vendor/benchmark")
set(SEQAN3_TEST_SRC_DIR "${PROJECT_BINARY_DIR}/vendor/googletest")

if(SEQAN3_SRC_DIR)
    find_package_message(SEQAN3_SRC_DIR "Found SEQAN3_SRC_DIR - ${SEQAN3_SRC_DIR}" "[${SEQAN3_SRC_DIR}]")
endif()

# required flags, includes, definitions and libraries for seqan3
set(SEQAN3_STRICT_CXX_FLAGS "-pedantic" "-Werror" "-Wall" "-Wextra")

# required flags, includes and libraries for seqan3/test/unit
set(SEQAN3_BENCHMARK_CXX_FLAGS ${SEQAN3_STRICT_CXX_FLAGS})
set(SEQAN3_BENCHMARK_INCLUDE_DIRS "")
set(SEQAN3_BENCHMARK_LIBRARIES "")
list(APPEND SEQAN3_BENCHMARK_INCLUDE_DIRS "${SEQAN3_BENCHMARK_SRC_DIR}/include/")

# required flags, includes and libraries for seqan3/test
set(SEQAN3_TEST_CXX_FLAGS ${SEQAN3_STRICT_CXX_FLAGS})
set(SEQAN3_TEST_INCLUDE_DIRS "")
set(SEQAN3_TEST_LIBRARIES "")
list(APPEND SEQAN3_TEST_INCLUDE_DIRS "${SEQAN3_TEST_SRC_DIR}/googletest/include/")

# commonly shared options:
set(SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "")
list(APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}")
list(APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")
list(APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}")
# force that libraries are installed to `lib/`, because GNUInstallDirs might install it into `lib64/`
list(APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_INSTALL_LIBDIR=${PROJECT_BINARY_DIR}/lib/")

macro(seqan3_require_ccache)
    find_program(CCACHE_PROGRAM ccache)
    find_package_message(CCACHE_PROGRAM_PRE "Finding program ccache" "[${CCACHE_PROGRAM}]")

    if(NOT CCACHE_PROGRAM)
        find_package_message(CCACHE_PROGRAM "Finding program ccache - Failed" "[${CCACHE_PROGRAM}]")
    else()
        find_package_message(CCACHE_PROGRAM "Finding program ccache - Success" "[${CCACHE_PROGRAM}]")
        if(CMAKE_VERSION VERSION_LESS 3.4)
            set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
            set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK "${CCACHE_PROGRAM}")
        else()
            # New option since cmake >= 3.4:
            # https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_LAUNCHER.html
            set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")

            # use ccache in external cmake projects
            list(APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_COMPILER_LAUNCHER=${CMAKE_CXX_COMPILER_LAUNCHER}")
        endif()
    endif()
    unset(CCACHE_PROGRAM)
endmacro()

macro (add_subdirectories)
    file (GLOB ENTRIES
          RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
          ${CMAKE_CURRENT_SOURCE_DIR}/[!.]*)

    foreach (ENTRY ${ENTRIES})
        if (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${ENTRY})
            if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${ENTRY}/CMakeLists.txt)
                add_subdirectory(${ENTRY})
            endif ()
        endif ()
    endforeach ()
    unset(ENTRIES)
endmacro ()
