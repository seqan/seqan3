#include <seqan3/std/ranges>
#include <vector>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/views/convert.hpp>
#include <seqan3/utility/views/to.hpp>

int main()
{
    // convert from int to bool
    std::vector<int> vec{7, 5, 0, 5, 0, 0, 4, 8, -3};

    // pipe notation
    seqan3::debug_stream << (vec | seqan3::views::convert<bool>) << '\n'; // [1,1,0,1,0,0,1,1,1]

    // combinability
    seqan3::debug_stream << (vec | seqan3::views::convert<bool> | std::views::reverse) << '\n'; // [1,1,1,0,0,1,0,1,1]

    // function notation and immediate conversion to vector again
    auto bool_vec = seqan3::views::convert<bool>(vec) | seqan3::views::to<std::vector<bool>>;
    seqan3::debug_stream << std::boolalpha << (bool_vec == std::vector<bool>{1, 1, 0, 1, 0, 0, 1, 1, 1})
                         << '\n'; // true
}
