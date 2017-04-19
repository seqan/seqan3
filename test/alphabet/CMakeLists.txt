add_subdirectory(aminoacid)
add_subdirectory(gap)
add_subdirectory(nucleotide)
add_subdirectory(quality)

add_executable(union_alphabet_test union_alphabet_test.cpp)
target_link_libraries(union_alphabet_test ${SEQAN3_LIBRARIES})
add_test(NAME alphabet_union_alphabet COMMAND union_alphabet_test)
