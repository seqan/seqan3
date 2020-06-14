# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

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
    add_subdirectories_of(${CMAKE_CURRENT_SOURCE_DIR})
endmacro ()
