#include <iostream>

#include <seqan3/range/view/take_until.hpp>
#include <seqan3/std/view/reverse.hpp>

using namespace seqan3;

int main()
{
//! [usage]
std::string vec{"foo\nbar"};
auto v = vec | view::take_until([] (char c) { return c == '\n'; });
std::cout << v << '\n'; // [f,o,o]

auto v2 = vec | view::reverse | view::take_until([] (char c) { return c == '\n'; });
std::cout << v2 << '\n'; // [r,a,b]
//! [usage]
}
