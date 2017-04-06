cmake_minimum_required (VERSION 3.0)
project (seqan3_alphabet CXX)

include_directories("../gaps")
include_directories("../nucleotide")
include_directories("../quality")

add_executable(test_illumina18 illumina18_test.cpp)

#add_dependencies(align_config_test gtest)
target_link_libraries(test_illumina18 ${SEQAN3_LIBRARIES})
add_test(NAME test_illumina18_test COMMAND test_illumina18)