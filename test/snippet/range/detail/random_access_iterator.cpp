#include <seqan3/range/detail/random_access_iterator.hpp>

namespace snippet
{
using namespace seqan3::detail;

//! [usage]
template <typename range_type>
class random_access_iterator : public random_access_iterator_base<range_type, random_access_iterator>
{
//...
};
//! [usage]
}

#if 0 // not supposed to compile
//! [not_usage]
template <typename range_type>
class random_access_iterator : public random_access_iterator_base<range_type, random_access_iterator<range_type>>
{
//...
};
//! [not_usage]
#endif

int main() {}
