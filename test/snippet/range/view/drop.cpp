#include <string>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/view/drop.hpp>               // provides view::drop
#include <seqan3/std/ranges>                        // provides std::view::reverse

int main()
{
    std::string s{"foobar"};
    auto v = s | seqan3::view::drop(3);
    seqan3::debug_stream << v << '\n';          // "bar"

    auto v2 = s | std::view::reverse | seqan3::view::drop(3);
    seqan3::debug_stream << v2 << '\n';         // "oof"
}
