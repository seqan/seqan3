#include <iostream>

#include <seqan3/io/stream/debug_stream.hpp>
#include <seqan3/range/view/take.hpp>
#include <seqan3/std/view/reverse.hpp>

using namespace seqan3;

int main()
{
//! [usage]
std::string vec{"foobar"};
auto v = vec | view::take(3);
debug_stream << v << '\n'; // [f,o,o]

auto v2 = vec | view::reverse | view::take(3);
debug_stream << v2 << '\n'; // [r,a,b]
//! [usage]
}
