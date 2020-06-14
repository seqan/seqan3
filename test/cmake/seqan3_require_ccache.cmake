# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

# Uses `ccache` to cache build results.
#
# See also
# * https://ccache.dev/
# * https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_LAUNCHER.html
macro (seqan3_require_ccache)
    find_program (CCACHE_PROGRAM ccache)
    find_package_message (CCACHE_PROGRAM_PRE "Finding program ccache" "[${CCACHE_PROGRAM}]")

    if (NOT CCACHE_PROGRAM)
        find_package_message (CCACHE_PROGRAM "Finding program ccache - Failed" "[${CCACHE_PROGRAM}]")
    else ()
        find_package_message (CCACHE_PROGRAM "Finding program ccache - Success" "[${CCACHE_PROGRAM}]")
        # New option since cmake >= 3.4:
        # https://cmake.org/cmake/help/latest/variable/CMAKE_LANG_COMPILER_LAUNCHER.html
        if (NOT CMAKE_VERSION VERSION_LESS 3.15) # cmake >= 3.15
            list (PREPEND CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
        else ()
            # prepend ccache to CMAKE_CXX_COMPILER_LAUNCHER
            list (INSERT CMAKE_CXX_COMPILER_LAUNCHER 0 "${CCACHE_PROGRAM}")
        endif ()

        # use ccache in external cmake projects
        list (APPEND SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS "-DCMAKE_CXX_COMPILER_LAUNCHER=${CMAKE_CXX_COMPILER_LAUNCHER}")
    endif ()
    unset (CCACHE_PROGRAM)
endmacro ()
