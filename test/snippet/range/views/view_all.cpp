#include <string>
#include <vector>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/views/view_all.hpp>       // provides views::all, attention not <seqan3/range/views/all.hpp>!

int main()
{
    std::string vec{"foobar"};
    auto v = vec | seqan3::views::all;   // pipe notation; v is of type std::string_view
    seqan3::debug_stream << v << '\n';   // "foobar"

    std::vector vec2{1, 2, 3};
    auto v2 = seqan3::views::all(vec2);  // functional notation; v is of type std::span
    seqan3::debug_stream << v2 << '\n';  // "[1, 2, 3]"
}
