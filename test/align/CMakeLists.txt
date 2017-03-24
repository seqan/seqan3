cmake_minimum_required(VERSION 3.2)
project(align_test CXX)

add_executable(alignment_class_test alignment_class_test.cpp)
add_dependencies(alignment_class_test gtest)
target_link_libraries(alignment_class_test ${SEQAN3_LIBRARIES})
add_test(NAME test_alignment_class_test COMMAND alignment_class_test)
