#include <seqan3/core/type_list.hpp>
#include <seqan3/test/testing_list.hpp>

using namespace seqan3;

template <typename T>
class my_typed_test : public ::testing::Test
{};

using types = type_list<int, char, double>;

TYPED_TEST_CASE(my_typed_test, as_testing_list<types>);
// same as TYPED_TEST_CASE(my_typed_test, ::testing::Types<int, char, double>);
