# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

include (diagnostics/list_missing_unit_tests)

# get_include_target (<VAR> TARGET dna4_test)
# get_include_target (<VAR> SOURCE "[test/unit/]alphabet/nucleotide/dna4_test.cpp")
# get_include_target (<VAR> HEADER "[include/seqan3/]alphabet/nucleotide/dna4.hpp")
function (get_include_target VAR source_type source)
    if (source_type STREQUAL "TARGET")
        # e.g. dna4_test.cpp
        get_target_property (target_source "${source}" SOURCES)
        # e.g. <seqan3-root>/test/unit/alphabet/nucleotide
        get_target_property (target_source_dir "${source}" SOURCE_DIR)
        # e.g. alphabet/nucleotide/dna4_test.cpp
        file (RELATIVE_PATH target_source_file "${CMAKE_SOURCE_DIR}" "${target_source_dir}/${target_source}")
    elseif (source_type STREQUAL "SOURCE")
        # e.g. alphabet/nucleotide/dna4_test.cpp
        set (target_source_file "${source}")
    elseif (source_type STREQUAL "HEADER")
        # e.g. alphabet/nucleotide
        get_filename_component (directory_name "${source}" DIRECTORY)
        # e.g. alphabet/nucleotide/dna4.hpp
        get_filename_component (filename_without_extension "${source}" NAME_WE)
        # e.g. alphabet/nucleotide/dna4_test.cpp
        set (target_source_file "${directory_name}/${filename_without_extension}_test.cpp")
    else ()
        message (STATUS "get_include_target: source_type ${source_type} unknown")
    endif ()

    string (REGEX REPLACE "\_test.cpp$" ".hpp" target_header_file "include/seqan3/${target_source_file}")
    seqan3_test_component (include_target "${target_header_file}" TARGET_UNIQUE_NAME)

    # e.g. include-seqan3-alphabet-nucleotide-dna4.hpp
    set (${VAR}
         "${include_target}.hpp"
         PARENT_SCOPE)
endfunction ()

function (add_include_target include_target)
    if (NOT TARGET ${include_target})
        # add_library(${include_target} STATIC IMPORTED)
        add_custom_target (${include_target})
        collect_include_target (${include_target})
    endif ()
endfunction ()

function (add_include_dependencies target target_cyclic_depending_includes)
    if (NOT SEQAN3_USE_INCLUDE_DEPENDENCIES)
        return ()
    endif ()

    if (NOT CMAKE_GENERATOR STREQUAL "Unix Makefiles")
        message (STATUS "add_include_dependencies: only works for Unix Makefiles")
        return ()
    endif ()

    if (NOT (DEFINED CMAKE_DEPENDS_USE_COMPILER) OR CMAKE_DEPENDS_USE_COMPILER)
        message (FATAL_ERROR "Starting with CMake 3.20, you need to specify -DCMAKE_DEPENDS_USE_COMPILER=OFF when "
                             "using -DSEQAN3_USE_INCLUDE_DEPENDENCIES=ON.")
    endif ()

    get_include_target (include_target TARGET "${target}")
    add_include_target ("${include_target}")
    add_dependencies (${include_target} ${target})

    # an include target can't depend on itself, we get a dependency cycle otherwise
    list (APPEND target_cyclic_depending_includes "${include_target}")

    # e.g. alphabet/nucleotide/dna4_test
    set (target_file "${CMAKE_CURRENT_BINARY_DIR}/${target}")
    # e.g. alphabet/nucleotide/dna4_test_dependencies_include.cmake
    set (target_dependencies_include_file "${target_file}_dependencies_include.cmake")
    # e.g. alphabet/nucleotide/dna4_test_dependencies.cmake
    set (target_dependencies_file "${target_file}_dependencies.cmake")

    file (WRITE "${target_dependencies_include_file}" #
          "if (EXISTS \"${target_dependencies_file}\")\n" #
          "  include (\"${target_dependencies_file}\")\n" #
          "endif ()\n")

    # e.g. alphabet/nucleotide/CMakeFiles/dna4_test.dir which is internally used by CMake Makefiles
    get_filename_component (target_internal_dir
                            "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/${target}.dir/${CMAKE_CFG_INTDIR}"
                            REALPATH)
    set (target_internal_dependency_info "${target_internal_dir}/DependInfo.cmake")
    set (target_internal_dependency_file "${target_internal_dir}/depend.internal")
    set (target_internal_dependency_make_file "${target_internal_dir}/depend.make")

    add_custom_command (OUTPUT "${target_internal_dependency_file}" "${target_dependencies_file}"
                        # generate alphabet/nucleotide/CMakeFiles/dna4_test.dir/depend.internal
                        COMMAND "${CMAKE_COMMAND}" #
                                "-E" #
                                "cmake_depends" #
                                "${CMAKE_GENERATOR}" #
                                "${CMAKE_SOURCE_DIR}" #
                                "${CMAKE_CURRENT_SOURCE_DIR}" #
                                "${CMAKE_BINARY_DIR}" #
                                "${CMAKE_CURRENT_BINARY_DIR}" #
                                "${target_internal_dependency_info}" #
                                "-DCMAKE_DEPENDS_USE_COMPILER=OFF"
                        # generate alphabet/nucleotide/dna4_test_dependencies.cmake
                        COMMAND "${CMAKE_COMMAND}" #
                                "-DTARGET=${target}" #
                                "-DTARGET_INTERNAL_DEPENDENCY_MAKE_FILE=${target_internal_dependency_make_file}" #
                                "-DSEQAN3_INCLUDE_DIR=${SEQAN3_INCLUDE_DIR}" #
                                "-DSEQAN3_TEST_CMAKE_MODULE_DIR=${SEQAN3_TEST_CMAKE_MODULE_DIR}" #
                                "-DTARGET_DEPENDENCIES_FILE=${target_dependencies_file}" #
                                "-DTARGET_CYCLIC_DEPENDING_INCLUDES=${target_cyclic_depending_includes}" #
                                "-P" #
                                "${SEQAN3_TEST_CMAKE_MODULE_DIR}/include_dependencies/generate_include_dependencies.cmake"
                        # touch alphabet/nucleotide/dna4_test_dependencies_include.cmake to "refresh" dependencies
                        COMMAND "${CMAKE_COMMAND}" #
                                "-E" #
                                "touch_nocreate" #
                                "${target_dependencies_include_file}"
                        BYPRODUCTS "${target_internal_dependency_file}" "${target_dependencies_file}"
                        DEPENDS "${SEQAN3_TEST_CMAKE_MODULE_DIR}/include_dependencies/generate_include_dependencies.cmake"
                        VERBATIM)

    add_custom_target ("${target}_dependencies" DEPENDS "${target_internal_dependency_file}")
    add_dependencies (${target} "${target}_dependencies")

    if (NOT TARGET all_dependencies)
        # custom target to build / collect all dependencies
        add_custom_target (all_dependencies)
    endif ()
    add_dependencies (all_dependencies "${target}_dependencies")

    include (${target_dependencies_include_file})
endfunction ()
