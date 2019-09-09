#include <string>
#include <vector>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/views/to_upper.hpp>

int main()
{
    std::string s{"CHanGED!"};
    std::vector<std::string> sv{"changed", "UNCHANGED!"};
    auto v1 = s | seqan3::views::to_upper;
    auto v2 = sv | seqan3::views::to_upper;

    seqan3::debug_stream << v1 << '\n'; // => "CHANGED!"
    seqan3::debug_stream << v2 << '\n'; // => ["CHANGED", "UNCHANGED!"]
}
