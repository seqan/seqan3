cmake_minimum_required(VERSION 3.2)
project(nucleotide_test CXX)

add_executable(dna5_test dna5_test.cpp)
add_executable(dna5_sequence_test dna5_sequence_test.cpp)
add_test(NAME test_dna5_test COMMAND dna5_test)
add_test(NAME test_dna5_sequence_test COMMAND dna5_sequence_test)

target_link_libraries(dna5_test ${SEQAN3_LIBRARIES})
target_link_libraries(dna5_sequence_test ${SEQAN3_LIBRARIES})
