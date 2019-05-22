#include <string>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/view/drop.hpp>               // provides view::drop
#include <seqan3/std/ranges>                        // provides std::view::reverse

using namespace seqan3;

int main()
{
    std::string s{"foobar"};
    auto v = s | view::drop(3);
    debug_stream << v << '\n';          // "bar"

    auto v2 = s | std::view::reverse | view::drop(3);
    debug_stream << v2 << '\n';         // "oof"
}
