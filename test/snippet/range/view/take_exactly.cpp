#include <string>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/view/take_exactly.hpp>           // provides view::take_exactly and view::take_exactly_or_throw

using namespace seqan3;

int main()
{
    std::string vec{"foobar"};
    auto v = vec | view::take_exactly(3);               // or view::take_exactly_or_throw
    debug_stream << v << '\n';                          // "foo"
    debug_stream << std::ranges::size(v) << '\n';       // 3


    auto v2 = vec | view::take_exactly(9);
    debug_stream << std::ranges::size(v2) << '\n';      // 9 <- here be dragons!
//     debug_stream << v2 << '\n';                      // possibly memory access violation

//     auto v3 = vec | view::take_exactly_or_throw(9);  // safe, throws immediately on construction

}
