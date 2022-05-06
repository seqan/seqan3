# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.10)

include (${CMAKE_SOURCE_DIR}/../cmake/seqan3_test_files.cmake)

function (copy_files_and_replace_seqan3 INPUT_PATH WILDCARDS OUTPUT_DIR_NAME NAME_TO_REPLACE_WITH_SEQAN3)

    # collect files to copy
    seqan3_test_files (input_files "${INPUT_PATH}" "${WILDCARDS}")

    # create output directory
    file (MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/additional_doc/${OUTPUT_DIR_NAME})

    # for each file: (1) read the file, (2) replace ${NAME_TO_REPLACE_WITH_SEQAN3} with "seqan3", (3) write the file
    foreach (input_filename ${input_files})
        file (READ "${INPUT_PATH}/${input_filename}" INPUT_CONTENT)
        string (REGEX REPLACE "${NAME_TO_REPLACE_WITH_SEQAN3}" "seqan3" OUTPUT_CONTENT ${INPUT_CONTENT})
        file (WRITE "${CMAKE_BINARY_DIR}/additional_doc/${OUTPUT_DIR_NAME}/${input_filename}" ${OUTPUT_CONTENT})
    endforeach ()

endfunction()

function (copy_tutorial_and_replace_seqan3 INPUT_PATH OUTPUT_DIR_NAME NAME_TO_REPLACE_WITH_SEQAN3)

    # collect files to copy
    seqan3_test_files (doc_index_md_files "${INPUT_PATH}" "index.md")

    # create output directory
    file (MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/additional_doc/${OUTPUT_DIR_NAME})

    # for each file: (1) read the file, (2) replace ${NAME_TO_REPLACE_WITH_SEQAN3} with "seqan3", (3) write the file
    foreach (index_md_filename ${doc_index_md_files})

        file (READ "${INPUT_PATH}/${index_md_filename}" INPUT_CONTENT)

        string (REGEX REPLACE "${NAME_TO_REPLACE_WITH_SEQAN3}" "seqan3" OUTPUT_CONTENT ${INPUT_CONTENT})

        # in case this is a documentation page, also change inlcudes and references.
        # since the regex should not match any content in snippets, it's fine like this.
        string(REGEX MATCHALL "[^ ]+\\.cpp" INCLUDED_CPP_FILES ${INPUT_CONTENT})

        foreach (cpp_file_path ${INCLUDED_CPP_FILES})
            get_filename_component(cpp_filename ${cpp_file_path} NAME)
            file (READ "${INPUT_PATH}../../${cpp_file_path}" CPP_CONTENT)
            string (REGEX REPLACE "${NAME_TO_REPLACE_WITH_SEQAN3}" "seqan3" CPP_OUTPUT_CONTENT ${CPP_CONTENT})
            file (WRITE "${CMAKE_BINARY_DIR}/additional_doc/${OUTPUT_DIR_NAME}/${cpp_filename}" ${CPP_OUTPUT_CONTENT})
        endforeach()

        string (REGEX REPLACE
                " [^ ]*/\([^/]*\)\.cpp"
                " ${CMAKE_BINARY_DIR}/additional_doc/${OUTPUT_DIR_NAME}/\\1.cpp"
                OUTPUT_CONTENT2
                ${OUTPUT_CONTENT})
        string (REGEX REPLACE "{#\(.*\)}" "{#seqan3_\\1}" OUTPUT_CONTENT3 ${OUTPUT_CONTENT2})

        file (WRITE "${CMAKE_BINARY_DIR}/additional_doc/index_files/${index_md_filename}" ${OUTPUT_CONTENT3})
    endforeach ()

endfunction()
