add_executable(char_test char_test.cpp)
target_link_libraries(char_test ${SEQAN3_LIBRARIES})
add_test(NAME alphabet_adaptation_char COMMAND char_test)

add_executable(uint_test uint_test.cpp)
target_link_libraries(uint_test ${SEQAN3_LIBRARIES})
add_test(NAME alphabet_adaptation_uint COMMAND uint_test)
