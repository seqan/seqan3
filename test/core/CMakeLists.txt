add_executable(pod_tuple_test pod_tuple_test.cpp)
target_link_libraries(pod_tuple_test ${SEQAN3_LIBRARIES})
add_test(NAME core_pod_tuple COMMAND pod_tuple_test)
