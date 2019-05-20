#include <iostream>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/view/take_line.hpp>
#include <seqan3/range/view/single_pass_input.hpp>
#include <seqan3/std/ranges>

using namespace seqan3;

int main()
{
{
//! [adaptor_def]
ranges::view::take_while([] (auto const & l) { return (l != '\r') && (l != '\n'); });
//! [adaptor_def]
}

{
//! [behaviour]
std::string vec{"foo\nbar"};
auto v = vec | view::take_line;
debug_stream << v << '\n'; // [f,o,o]

auto v2 = vec | std::view::reverse | view::take_line;
debug_stream << v2 << '\n'; // [r,a,b]
debug_stream << v2 << '\n'; // [r,a,b] (parsing it again gives us the same result)
//! [behaviour]
}

{
//! [tokenise]
std::string vec{"foo\nbar"};
auto v = vec | view::single_pass_input | view::take_line;
debug_stream << v << '\n'; // [f,o,o]
debug_stream << v << '\n'; // [b,a,r] (parsing it again gives us the next line)
//! [tokenise]
}
}
