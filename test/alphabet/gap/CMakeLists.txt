add_executable(gap_test gap_test.cpp)
add_test(NAME test_gap_test COMMAND gap_test)
target_link_libraries(gap_test ${SEQAN3_LIBRARIES})

add_executable(gapped_alphabet_test gapped_alphabet_test.cpp)
add_test(NAME test_gapped_alphabet_test COMMAND gapped_alphabet_test)
target_link_libraries(gapped_alphabet_test ${SEQAN3_LIBRARIES})
