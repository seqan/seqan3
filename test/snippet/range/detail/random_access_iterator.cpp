#include <seqan3/core/range/detail/random_access_iterator.hpp>

template <typename range_type>
class my_random_access_iterator : public seqan3::detail::random_access_iterator_base<range_type, seqan3::detail::random_access_iterator>
{
//...
};

int main() {}
