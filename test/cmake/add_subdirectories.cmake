# SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cmake_minimum_required (VERSION 3.10)

# Calls add_subdirectory on all (direct) subdirectories of the given directory if they contain a `CMakeLists.txt`
#
# Example:
# If we have
# * /some/path/directory/subdir1/CMakeLists.txt
# * /some/path/directory/subdir2/CMakeLists.txt
# * /some/path/directory/subdir3/CMakeLists.txt
# * /some/path/directory/subdir4/has-no-CMakeLists.txt
#
# This macro calls
# * `add_subdirectory ("/some/path/directory/subdir1")`,
# * `add_subdirectory ("/some/path/directory/subdir2")`, and
# * `add_subdirectory ("/some/path/directory/subdir3")`,
# but not `add_subdirectory ("/some/path/directory/subdir4")`, because it does not contain a CMakeLists.txt
macro (add_subdirectories_of directory)
    file (GLOB ENTRIES
          RELATIVE ${directory}
          ${directory}/[!.]*)

    foreach (ENTRY ${ENTRIES})
        if (IS_DIRECTORY ${directory}/${ENTRY})
            if (EXISTS ${directory}/${ENTRY}/CMakeLists.txt)
                add_subdirectory (${directory}/${ENTRY} ${CMAKE_CURRENT_BINARY_DIR}/${ENTRY})
            endif ()
        endif ()
    endforeach ()
    unset (ENTRIES)
endmacro ()

macro (add_subdirectories)
    add_subdirectories_of (${CMAKE_CURRENT_SOURCE_DIR})
endmacro ()
