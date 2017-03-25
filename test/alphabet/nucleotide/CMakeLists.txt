cmake_minimum_required(VERSION 3.2)
project(nucleotide_test CXX)

add_executable(dna5_test dna5_test.cpp)
add_test(NAME test_dna5_test COMMAND dna5_test)
target_link_libraries(dna5_test ${SEQAN3_LIBRARIES})

add_executable(dna5_sequence_test dna5_sequence_test.cpp)
add_test(NAME test_dna5_sequence_test COMMAND dna5_sequence_test)
target_link_libraries(dna5_sequence_test ${SEQAN3_LIBRARIES})

add_executable(alphabet_nucleotide_nucl16_test nucl16_test.cpp)
add_test(NAME test_alphabet_nucleotide_nucl16_test COMMAND alphabet_nucleotide_nucl16_test)
target_link_libraries(alphabet_nucleotide_nucl16_test ${SEQAN3_LIBRARIES})

add_executable(alphabet_nucleotide_nucl16_sequence_test nucl16_sequence_test.cpp)
add_test(NAME test_alphabet_nucleotide_nucl16_sequence_test COMMAND alphabet_nucleotide_nucl16_sequence_test)
target_link_libraries(alphabet_nucleotide_nucl16_sequence_test ${SEQAN3_LIBRARIES})

add_executable(rna4_test rna4_test.cpp)
add_test(NAME test_rna4_test COMMAND rna4_test)
target_link_libraries(rna4_test ${SEQAN3_LIBRARIES})

add_executable(rna4_container_test rna4_container_test.cpp)
add_test(NAME test_rna4_container_test COMMAND rna4_container_test)
target_link_libraries(rna4_container_test ${SEQAN3_LIBRARIES})

