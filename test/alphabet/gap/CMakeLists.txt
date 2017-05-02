add_executable(gap_test gap_test.cpp)
target_link_libraries(gap_test ${SEQAN3_LIBRARIES})
add_test(NAME alphabet_gap_gap COMMAND gap_test)

add_executable(gapped_alphabet_test gapped_alphabet_test.cpp)
target_link_libraries(gapped_alphabet_test ${SEQAN3_LIBRARIES})
add_test(NAME alphabet_gap_gapped_alphabet COMMAND gapped_alphabet_test)
