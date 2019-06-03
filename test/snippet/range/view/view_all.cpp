#include <string>
#include <vector>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/view/view_all.hpp>       // provides view::all, attention not <seqan3/range/view/all.hpp>!

using namespace seqan3;

int main()
{
    std::string vec{"foobar"};
    auto v = vec | view::all;   // pipe notation; v is of type std::string_view
    debug_stream << v << '\n';  // "foobar"

    std::vector vec2{1, 2, 3};
    auto v2 = view::all(vec2);  // functional notation; v is of type std::span
    debug_stream << v2 << '\n'; // "[1, 2, 3]"
}
