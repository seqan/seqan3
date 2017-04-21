add_executable(int_types_test int_types_test.cpp)
target_link_libraries(int_types_test ${SEQAN3_LIBRARIES})
add_test(NAME core_detail_int_types COMMAND int_types_test)
