#include <iostream>

#include <seqan3/range/view/take_exactly.hpp>

using namespace seqan3;

int main()
{

{
//! [usage]
std::string vec{"foobar"};
auto v = vec | view::take_exactly(3);        // or view::take_exactly_or_throw
std::cout << v << '\n';                      // [f,o,o]
std::cout << ranges::size(v) << v << '\n';   // 3
//! [usage]
}

{
//! [shorter_sequence]
std::string vec{"foo"};
auto v = vec | view::take_exactly(4);
std::cout << v << '\n';                          // [f,o,o]
std::cout << ranges::size(v) << v << '\n';       // 4 <- here be dragons!

auto v2 = vec | view::take_exactly_or_throw(4);  // throws immediately on construction
//! [shorter_sequence]
}
}
