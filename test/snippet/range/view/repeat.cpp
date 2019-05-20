#include <iostream>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/view/repeat.hpp>
#include <seqan3/std/ranges>

using namespace seqan3;

int main()
{
    auto v = view::repeat('A');

    debug_stream << *std::ranges::begin(v) << std::endl; // prints 'A'
    debug_stream << v[12355] << std::endl;               // also prints 'A'. It always prints 'A'

    v[1345] = 'C';

    // Now it always prints 'C'
    debug_stream << *std::ranges::begin(v) << std::endl; // prints 'C'
    debug_stream << v[12355] << std::endl;               // prints 'C'
}
