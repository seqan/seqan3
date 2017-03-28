cmake_minimum_required(VERSION 3.0)
project(container_tests CXX)

add_executable(ra_iterator_test ra_iterator_test.cpp)
add_test(NAME test_ra_iterator_test COMMAND ra_iterator_test)
target_link_libraries(ra_iterator_test ${SEQAN3_LIBRARIES})
