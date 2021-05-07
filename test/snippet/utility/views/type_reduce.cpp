#include <string>
#include <vector>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/views/type_reduce.hpp>

int main()
{
    std::string const vec{"foobar"};
    auto v = vec | seqan3::views::type_reduce;          // pipe notation; v is of type std::string_view
    seqan3::debug_stream << v << '\n';                  // "foobar"

    std::vector vec2{1, 2, 3};
    auto v2 = seqan3::views::type_reduce(vec2);         // functional notation; v2 is of type std::span
    seqan3::debug_stream << v2 << '\n';                 // "[1, 2, 3]"
}
