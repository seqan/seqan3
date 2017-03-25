cmake_minimum_required (VERSION 3.0)
project (seqan3_alphabet CXX)

include_directories("../gaps")
include_directories("../nucleotide")
include_directories("../quality")

add_executable(test_phred phred_test.cpp)

#add_dependencies(align_config_test gtest)
target_link_libraries(test_phred ${SEQAN3_LIBRARIES})
add_test(NAME test_phred_test COMMAND test_phred)