cmake_minimum_required (VERSION 3.0)
project (seqan3_alphabet CXX)

add_subdirectory(detail)

add_executable(test_pod_tuple pod_tuple_test.cpp)

target_link_libraries(test_pod_tuple ${SEQAN3_LIBRARIES})
add_test(NAME test_pod_tuple_test COMMAND test_pod_tuple)
