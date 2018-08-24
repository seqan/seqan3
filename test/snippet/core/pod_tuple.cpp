#include <seqan3/core/pod_tuple.hpp>
#include <iostream>

using namespace seqan3;

int main()
{
//! [usage]
pod_tuple<int, float> t{3, 4.7};
static_assert(std::is_pod_v<pod_tuple<int, float>>);

// template parameters are automatically deduced:
pod_tuple t2{17, 3.7f, 19l};

std::cout << std::get<0>(t2) << '\n'; // 17

auto [ i, f, l ] = t2; // creates an int i with value 17, float f...
//! [usage]
(void) t;
(void) i;
(void) f;
(void) l;
}
