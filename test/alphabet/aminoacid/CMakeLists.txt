add_executable(aa27_test aa27_test.cpp)
target_link_libraries(aa27_test ${SEQAN3_LIBRARIES})
add_test(NAME alphabet_aminoacid_aa27 COMMAND aa27_test)

add_executable(aa27_container_test aa27_container_test.cpp)
target_link_libraries(aa27_container_test ${SEQAN3_LIBRARIES})
add_test(NAME alphabet_aminoacid_aa27_container COMMAND aa27_container_test)
