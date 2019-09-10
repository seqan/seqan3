#include <string>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/views/slice.hpp>     // provides views::slice
#include <seqan3/std/ranges>                // provides std::views::reverse

int main()
{
    std::string s{"foobar"};
    auto v = s | seqan3::views::slice(1,4);
    seqan3::debug_stream << v << '\n';      // "oob"

    auto v2 = s | std::views::reverse | seqan3::views::slice(1, 4);
    seqan3::debug_stream << v2 << '\n';     // "abo"

}
