#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/views/repeat.hpp>
#include <seqan3/range/views/zip.hpp>

int main()
{
    std::vector<size_t> vec1{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    std::vector<std::string> vec2{"hello", "goodbye", "hello", "goodbye"};
    auto rep = seqan3::views::repeat('L');

    auto v = seqan3::views::zip(vec1, vec2, rep);
    seqan3::debug_stream << v << '\n';  // Prints [(0,hello,L),(1,goodbye,L),(2,hello,L),(3,goodbye,L)]

    *v.begin() = std::tuple{9, "MOO", 'P'};
    seqan3::debug_stream << v << '\n';  // Prints [(9,MOO,P),(1,goodbye,P),(2,hello,P),(3,goodbye,P)]

    auto v2 = v | std::views::reverse;
    seqan3::debug_stream << v2 << '\n'; // Prints [(3,goodbye,P),(2,hello,P),(1,goodbye,P),(9,MOO,P)]
}
