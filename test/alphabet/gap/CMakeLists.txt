add_executable(gap_test gap_test.cpp)
add_test(NAME test_gap_test COMMAND gap_test)
target_link_libraries(gap_test ${SEQAN3_LIBRARIES})
