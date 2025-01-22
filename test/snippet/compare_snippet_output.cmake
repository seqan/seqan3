# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

include ("${CMAKE_CURRENT_LIST_DIR}/../cmake/seqan3_path_longest_stem.cmake")

message (STATUS "TARGET_FILE: ${TARGET_FILE}")
message (STATUS "SOURCE_FILE: ${SOURCE_FILE}")

get_filename_component (TARGET_FILE_DIR "${TARGET_FILE}" DIRECTORY)
seqan3_path_longest_stem (SOURCE_FILE_NAME "${SOURCE_FILE}")
get_filename_component (SOURCE_FILE_DIR "${SOURCE_FILE}" DIRECTORY)

set (ACTUAL_OUTPUT_FILE "${TARGET_FILE}.out")
set (EXPECTED_OUTPUT_FILE "${SOURCE_FILE_DIR}/${SOURCE_FILE_NAME}.out")
set (ACTUAL_ERROR_FILE "${TARGET_FILE}.err")
set (EXPECTED_ERROR_FILE "${SOURCE_FILE_DIR}/${SOURCE_FILE_NAME}.err")

message (STATUS "ACTUAL_OUTPUT_FILE: ${ACTUAL_OUTPUT_FILE}")
message (STATUS "EXPECTED_OUTPUT_FILE: ${EXPECTED_OUTPUT_FILE}")

message (STATUS "ACTUAL_ERROR_FILE: ${ACTUAL_ERROR_FILE}")
message (STATUS "EXPECTED_ERROR_FILE: ${EXPECTED_ERROR_FILE}")

# execute snippet
execute_process (COMMAND "${TARGET_FILE}"
                 WORKING_DIRECTORY "${TARGET_FILE_DIR}"
                 OUTPUT_FILE "${ACTUAL_OUTPUT_FILE}"
                 ERROR_FILE "${ACTUAL_ERROR_FILE}"
                 RESULT_VARIABLE error_result)

if (error_result) # != 0 return code
    message (SEND_ERROR "error: executing snippet exited with '${error_result}'")
endif ()

function (compare_files actual_file expected_file)
    file (READ "${actual_file}" actual_output)

    if (actual_output AND EXISTS "${expected_file}")
        execute_process (COMMAND ${CMAKE_COMMAND} -E compare_files "${actual_file}" "${expected_file}"
                         RESULT_VARIABLE error_result)

        if (NOT error_result) # == 0 return code => files are identical
            # if successful move one
            return ()
        endif ()

        message (SEND_ERROR "error: `${actual_file}` and `${expected_file}` differ (exited with '${error_result}')")

        find_package (Git)
        if (Git_FOUND)
            execute_process (COMMAND "${GIT_EXECUTABLE}" diff --no-index "${expected_file}" "${actual_file}")
        endif ()
    elseif (EXISTS "${expected_file}")
        message (SEND_ERROR "error: `${expected_file}` exists, but `${actual_file}` has no output.")
    elseif (actual_output)
        message (SEND_ERROR "error: `${actual_file}` has output, but `${expected_file}` does not exist.")
    else ()
        message (STATUS "Output matches ${expected_file}")
    endif ()
endfunction ()

# compare output to file
compare_files ("${ACTUAL_OUTPUT_FILE}" "${EXPECTED_OUTPUT_FILE}")
compare_files ("${ACTUAL_ERROR_FILE}" "${EXPECTED_ERROR_FILE}")
