# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

# Finds all files relative to the `test_base_path_` which satisfy the given file pattern.
#
# Example:
# seqan3_test_files (header_files "/seqan3/include" "*.hpp;*.h")
#
# The variable `header_files` will contain:
#   seqan3/alphabet/adaptation/all.hpp
#   seqan3/alphabet/adaptation/char.hpp
#   seqan3/alphabet/adaptation/concept.hpp
#   seqan3/alphabet/adaptation/uint.hpp
#   seqan3/alphabet/all.hpp
#   seqan3/alphabet/dna5_detail.hpp
#   ....
macro (seqan3_test_files VAR test_base_path_ extension_wildcards)
    # test_base_path is /home/.../seqan3/test/
    get_filename_component (test_base_path "${test_base_path_}" ABSOLUTE)
    file (RELATIVE_PATH test_base_path_relative "${CMAKE_CURRENT_SOURCE_DIR}" "${test_base_path}")
    # ./ is a hack to deal with empty test_base_path_relative
    set (test_base_path_relative "./${test_base_path_relative}")
    # collect all cpp files
    set (${VAR} "")
    foreach (extension_wildcard ${extension_wildcards})
        file (GLOB_RECURSE test_files
              RELATIVE "${test_base_path}"
              "${test_base_path_relative}/${extension_wildcard}")
        list (APPEND ${VAR} ${test_files})
    endforeach ()

    unset (test_base_path)
    unset (test_base_path_relative)
endmacro ()
