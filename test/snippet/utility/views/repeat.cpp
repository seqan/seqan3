#include <seqan3/std/ranges>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/views/repeat.hpp>

int main()
{
    auto v = seqan3::views::repeat('A');

    seqan3::debug_stream << *std::ranges::begin(v) << '\n'; // prints 'A'
    seqan3::debug_stream << v[12355] << '\n';               // also prints 'A'. It always prints 'A'

    v[1345] = 'C';

    // Now it always prints 'C'
    seqan3::debug_stream << *std::ranges::begin(v) << '\n'; // prints 'C'
    seqan3::debug_stream << v[12355] << '\n';               // prints 'C'
}
