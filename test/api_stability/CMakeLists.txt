cmake_minimum_required(VERSION 3.8)

project (api_stability)

find_package(Git)

include (../seqan3-test.cmake)

include (ExternalProject)
set (SEQAN3_LAST_STABLE_VERSION "3.0.3")
set (SEQAN3_LAST_STABLE_SOURCE_URL "https://github.com/seqan/seqan3/releases/download/${SEQAN3_LAST_STABLE_VERSION}/seqan3-${SEQAN3_LAST_STABLE_VERSION}-Source.tar.xz")
set (SEQAN3_LAST_STABLE_SOURCE_SHA256 "42fca82e7e03ba2ed2026be739cd5e0be785152f872fdd0bac13b6f836035fc0")

ExternalProject_Add(api_stability_unit
  URL               ${SEQAN3_LAST_STABLE_SOURCE_URL}
  URL_HASH          SHA256=${SEQAN3_LAST_STABLE_SOURCE_SHA256}
  SOURCE_SUBDIR     "test/unit"
  CMAKE_ARGS        ${SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS}
                    -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -DSEQAN3_REMOVE_DEPRECATED_310=1"
                    -DCMAKE_CXX_COMPILER_LAUNCHER=${CMAKE_CXX_COMPILER_LAUNCHER}
                    -DCMAKE_PREFIX_PATH=${CMAKE_CURRENT_LIST_DIR}/../../build_system # use current build_system
  PATCH_COMMAND     ${CMAKE_COMMAND} -E echo "<SOURCE_DIR>"
  COMMAND           ${CMAKE_COMMAND} -E remove_directory "<SOURCE_DIR>/include"
  COMMAND           ${CMAKE_COMMAND} -E remove_directory "<SOURCE_DIR>/submodules"
  COMMAND           ${CMAKE_COMMAND} -E remove_directory "<SOURCE_DIR>/build_system"
  COMMAND           sh -c "${GIT_EXECUTABLE} apply --verbose ${CMAKE_CURRENT_LIST_DIR}/${SEQAN3_LAST_STABLE_VERSION}/*.patch"
  TEST_BEFORE_INSTALL TRUE
  INSTALL_COMMAND   ""
)

# Same as above, but with snippet instead of unit
ExternalProject_Add(api_stability_snippet
  URL               ${SEQAN3_LAST_STABLE_SOURCE_URL}
  URL_HASH          SHA256=${SEQAN3_LAST_STABLE_SOURCE_SHA256}
  SOURCE_SUBDIR     "test/snippet"
  CMAKE_ARGS        ${SEQAN3_EXTERNAL_PROJECT_CMAKE_ARGS}
                    -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -DSEQAN3_REMOVE_DEPRECATED_310=1"
                    -DCMAKE_CXX_COMPILER_LAUNCHER=${CMAKE_CXX_COMPILER_LAUNCHER}
                    -DCMAKE_PREFIX_PATH=${CMAKE_CURRENT_LIST_DIR}/../../build_system # use current build_system
  PATCH_COMMAND     ${CMAKE_COMMAND} -E echo "<SOURCE_DIR>"
  COMMAND           ${CMAKE_COMMAND} -E remove_directory "<SOURCE_DIR>/include"
  COMMAND           ${CMAKE_COMMAND} -E remove_directory "<SOURCE_DIR>/submodules"
  COMMAND           ${CMAKE_COMMAND} -E remove_directory "<SOURCE_DIR>/build_system"
  COMMAND           sh -c "${GIT_EXECUTABLE} apply -v ${CMAKE_CURRENT_LIST_DIR}/${SEQAN3_LAST_STABLE_VERSION}/*.patch"
  TEST_BEFORE_INSTALL TRUE
  INSTALL_COMMAND   ""
)
