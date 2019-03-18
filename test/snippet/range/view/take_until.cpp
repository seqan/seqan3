#include <iostream>

#include <seqan3/io/stream/debug_stream.hpp>
#include <seqan3/range/view/take_until.hpp>
#include <seqan3/std/ranges>

using namespace seqan3;

int main()
{
//! [usage]
std::string vec{"foo\nbar"};
auto v = vec | view::take_until([] (char c) { return c == '\n'; });
debug_stream << v << '\n'; // [f,o,o]

auto v2 = vec | std::view::reverse | view::take_until([] (char c) { return c == '\n'; });
debug_stream << v2 << '\n'; // [r,a,b]
//! [usage]
}
