#include <string>
#include <vector>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/views/to_lower.hpp>

int main()
{
    std::string s{"CHanGED!"};
    std::vector<std::string> sv{"CHANGED", "unchanged!"};
    auto v1 = s | seqan3::views::to_lower;
    auto v2 = sv | seqan3::views::to_lower;

    seqan3::debug_stream << v1 << '\n'; // => "changed!"
    seqan3::debug_stream << v2 << '\n'; // => ["changed", "unchanged!"]
}
