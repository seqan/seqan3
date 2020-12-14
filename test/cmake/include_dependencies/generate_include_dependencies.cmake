# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.7)

function (generate_include_dependencies_impl)
    cmake_parse_arguments(
        ""
        ""
        "TARGET;TARGET_INTERNAL_DEPENDENCY_MAKE_FILE;SEQAN3_INCLUDE_DIR;TARGET_DEPENDENCIES_FILE"
        "TARGET_CYCLIC_DEPENDING_INCLUDES"
        ${ARGN}
    )

    if (NOT EXISTS "${_TARGET_INTERNAL_DEPENDENCY_MAKE_FILE}")
        return ()
    endif ()

    # read in file and filter out linebreaks
    file (STRINGS "${_TARGET_INTERNAL_DEPENDENCY_MAKE_FILE}" header_files)
    list (FILTER header_files INCLUDE REGEX "${_SEQAN3_INCLUDE_DIR}/seqan3")
    if (CMAKE_VERSION VERSION_LESS 3.12)
        set (_header_files "${header_files}")
        set (header_files "")
        foreach (header_file ${_header_files})
            string(REGEX REPLACE "^.+: " "" header_file "${header_file}")
            list (APPEND header_files "${header_file}")
        endforeach ()
    else () # ^^^ workaround / no workaround vvv
        # list (TRANSFORM) needs cmake >= 3.12
        list (TRANSFORM header_files REPLACE "^.+: " "")
    endif ()

    set (script "")
    foreach (header_file ${header_files})
        file (RELATIVE_PATH relative_header_file "${_SEQAN3_INCLUDE_DIR}/seqan3" "${header_file}")
        get_include_target (include_depended_target HEADER "${relative_header_file}")

        string (APPEND script "add_include_target (\"${include_depended_target}\")\n")
        # exclude known cyclic header
        if (NOT "${include_depended_target}" IN_LIST _TARGET_CYCLIC_DEPENDING_INCLUDES)
            string (APPEND script "add_dependencies (${_TARGET} \"${include_depended_target}\")\n")
        else ()
            string (APPEND script "# add_dependencies (${_TARGET} \"${include_depended_target}\") # known cycle\n")
        endif ()
        string (APPEND script "\n")
    endforeach ()

    file (WRITE "${_TARGET_DEPENDENCIES_FILE}" "${script}")
endfunction ()

if (CMAKE_SCRIPT_MODE_FILE)
    list (APPEND CMAKE_MODULE_PATH "${SEQAN3_TEST_CMAKE_MODULE_DIR}")
    include (include_dependencies/add_include_dependencies)
    include (seqan3_test_component)

    message (STATUS "Generate include dependencies of target ${TARGET}")

    generate_include_dependencies_impl (
        # e.g. dna4_test
        TARGET "${TARGET}"
        # e.g. alphabet/nucleotide/CMakeFiles/dna4_test.dir/depend.make
        TARGET_INTERNAL_DEPENDENCY_MAKE_FILE "${TARGET_INTERNAL_DEPENDENCY_MAKE_FILE}"
        # e.g. /seqan3-repo/include
        SEQAN3_INCLUDE_DIR "${SEQAN3_INCLUDE_DIR}"
        # e.g. alphabet/nucleotide/dna4_test_dependencies.cmake (will be generated)
        TARGET_DEPENDENCIES_FILE "${TARGET_DEPENDENCIES_FILE}"
        TARGET_CYCLIC_DEPENDING_INCLUDES "${TARGET_CYCLIC_DEPENDING_INCLUDES}"
    )
endif ()
