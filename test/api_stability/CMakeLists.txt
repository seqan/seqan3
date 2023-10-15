cmake_minimum_required (VERSION 3.8)

project (api_stability)

find_package (Git)

include (../seqan3-test.cmake)

include (ExternalProject)
set (SEQAN3_LAST_STABLE_VERSION "3.3.0")
set (SEQAN3_LAST_STABLE_SOURCE_URL
     "https://github.com/seqan/seqan3/releases/download/${SEQAN3_LAST_STABLE_VERSION}/seqan3-${SEQAN3_LAST_STABLE_VERSION}-Source.tar.xz"
)
set (SEQAN3_LAST_STABLE_SOURCE_SHA256 "da2fb621268ebc52b9cc26087e96f4a94109db1f4f28d363d19c7c9cdbd788b1")

set (SEQAN3_API_STABILITY_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSEQAN3_DISABLE_DEPRECATED_WARNINGS=1")

macro (seqan3_api_stability_test test_suite git_apply)
    ExternalProject_Add (
        api_stability_${test_suite}
        URL ${SEQAN3_LAST_STABLE_SOURCE_URL}
        URL_HASH SHA256=${SEQAN3_LAST_STABLE_SOURCE_SHA256}
        SOURCE_SUBDIR "test/${test_suite}"
        CMAKE_ARGS ${SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS}
                   -DCMAKE_CXX_FLAGS=${SEQAN3_API_STABILITY_CXX_FLAGS}
                   -DCMAKE_CXX_COMPILER_LAUNCHER=${CMAKE_CXX_COMPILER_LAUNCHER}
                   -DCMAKE_PREFIX_PATH=${CMAKE_CURRENT_LIST_DIR}/../../build_system # use current build_system
                   -DSEQAN3_GENERATE_SNIPPETS=OFF
                   -DSEQAN3_VERBOSE_TESTS=OFF
        PATCH_COMMAND ${CMAKE_COMMAND} -E echo "<SOURCE_DIR>"
        COMMAND ${CMAKE_COMMAND} -E remove_directory "<SOURCE_DIR>/include"
        COMMAND ${CMAKE_COMMAND} -E remove_directory "<SOURCE_DIR>/submodules"
        COMMAND ${CMAKE_COMMAND} -E remove_directory "<SOURCE_DIR>/build_system"
        COMMAND sh -c "${git_apply}"
        TEST_BEFORE_INSTALL TRUE
        TEST_COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -j2
        INSTALL_COMMAND "")
endmacro ()

if (EXISTS "${CMAKE_CURRENT_LIST_DIR}/${SEQAN3_LAST_STABLE_VERSION}")
    set (SEQAN3_GIT_APPLY
         "${GIT_EXECUTABLE} apply --verbose ${CMAKE_CURRENT_LIST_DIR}/${SEQAN3_LAST_STABLE_VERSION}/*.patch")
else ()
    set (SEQAN3_GIT_APPLY "true")
    message (STATUS "No patches found in ${CMAKE_CURRENT_LIST_DIR}/${SEQAN3_LAST_STABLE_VERSION}.")
endif ()

seqan3_api_stability_test ("snippet" "${SEQAN3_GIT_APPLY}")
seqan3_api_stability_test ("unit" "${SEQAN3_GIT_APPLY}")
