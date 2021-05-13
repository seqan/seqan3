#include <seqan3/std/ranges>
#include <string>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/detail/take_view.hpp>

#ifdef SEQAN3_DEPRECATED_310
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
int main()
{
    std::string vec{"foobar"};
    auto v = vec | seqan3::views::take(3);
    seqan3::debug_stream << v << '\n'; // [f,o,o]

    auto v2 = vec | std::views::reverse | seqan3::views::take(3);
    seqan3::debug_stream << v2 << '\n'; // [r,a,b]
}
#pragma GCC diagnostic pop
#endif // SEQAN3_DEPRECATED_310
