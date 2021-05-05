#include <seqan3/std/ranges>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/utility/views/convert.hpp>
#include <seqan3/utility/views/to.hpp>

int main()
{
    // convert from int to bool
    std::vector<int>  vec{7, 5, 0, 5, 0, 0, 4, 8, -3};

    // pipe notation
    auto v = vec | seqan3::views::convert<bool>; // == [1, 1, 0, 1, 0, 0, 1, 1, 1];

    // function notation and immediate conversion to vector again
    auto v2 = seqan3::views::convert<bool>(vec) | seqan3::views::to<std::vector<bool>>;

    // combinability
    auto v3 = vec | seqan3::views::convert<bool> | std::views::reverse; // == [1, 1, 1, 0, 0, 1, 0, 1, 1];
}
