#include <seqan3/std/ranges>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/utility/views/deep.hpp>

namespace my
{
// You can create a permanent alias:
inline auto const deep_reverse = seqan3::views::deep{std::views::reverse};
}

int main()
{
    using namespace seqan3::literals;

    std::vector<seqan3::dna5_vector> foo{"AAATTT"_dna5, "CCCGGG"_dna5};

    auto r = foo | std::views::reverse;                     // == [ [C,C,C,G,G,G], [A,A,A,T,T,T] ]

    // These two are equivalent:
    auto e = foo | my::deep_reverse;                       // == [ [T,T,T,A,A,A], [G,G,G,C,C,C] ]
    auto d = foo | seqan3::views::deep{std::views::reverse}; // == [ [T,T,T,A,A,A], [G,G,G,C,C,C] ]
}
