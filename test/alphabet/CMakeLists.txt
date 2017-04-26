cmake_minimum_required(VERSION 3.2)
project(alphabet_test CXX)

add_all_subdirectories()

add_executable(alphabet_union_alphabet_test union_alphabet.cpp)
add_test(NAME test_alphabet_union_alphabet_test COMMAND alphabet_union_alphabet_test)
target_link_libraries(alphabet_union_alphabet_test ${SEQAN3_LIBRARIES})
