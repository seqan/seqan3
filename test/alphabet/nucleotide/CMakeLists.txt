add_executable(dna4_test dna4_test.cpp)
target_link_libraries(dna4_test ${SEQAN3_LIBRARIES})
add_test(NAME alphabet_nucleotide_dna4 COMMAND dna4_test)

add_executable(rna4_test rna4_test.cpp)
target_link_libraries(rna4_test ${SEQAN3_LIBRARIES})
add_test(NAME alphabet_nucleotide_rna4 COMMAND rna4_test)

add_executable(dna5_test dna5_test.cpp)
target_link_libraries(dna5_test ${SEQAN3_LIBRARIES})
add_test(NAME alphabet_nucleotide_dna5 COMMAND dna5_test)

add_executable(rna5_test rna5_test.cpp)
target_link_libraries(rna5_test ${SEQAN3_LIBRARIES})
add_test(NAME alphabet_nucleotide_rna5 COMMAND rna5_test)

add_executable(nucl16_test nucl16_test.cpp)
target_link_libraries(nucl16_test ${SEQAN3_LIBRARIES})
add_test(NAME alphabet_nucleotide_nucl16 COMMAND nucl16_test)
