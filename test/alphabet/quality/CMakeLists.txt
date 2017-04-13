cmake_minimum_required (VERSION 3.0)
project (seqan3_alphabet CXX)

include_directories("../gaps")
include_directories("../nucleotide")
include_directories("../quality")

add_executable(test_illumina18 illumina18_test.cpp)
target_link_libraries(test_illumina18 ${SEQAN3_LIBRARIES})
add_test(NAME test_illumina18_test COMMAND test_illumina18)

add_executable(test_quality_composition composition_test.cpp)
target_link_libraries(test_quality_composition ${SEQAN3_LIBRARIES})
add_test(NAME test_quality_composition_test COMMAND test_quality_composition)
