#include <string>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/view/take_exactly.hpp>           // provides view::take_exactly and view::take_exactly_or_throw

int main()
{
    std::string vec{"foobar"};
    auto v = vec | seqan3::view::take_exactly(3);               // or seqan3::view::take_exactly_or_throw
    seqan3::debug_stream << v << '\n';                          // "foo"
    seqan3::debug_stream << std::ranges::size(v) << '\n';       // 3


    auto v2 = vec | seqan3::view::take_exactly(9);
    seqan3::debug_stream << std::ranges::size(v2) << '\n';      // 9 <- here be dragons!

}
