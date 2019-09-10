#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/range/views/deep.hpp>
#include <seqan3/std/ranges>

namespace my
{
inline auto const deep_take1 = seqan3::views::deep{std::views::take(1)};
}

int main()
{
    using seqan3::operator""_dna5;

    std::vector<seqan3::dna5_vector> foo{"AAATTT"_dna5, "CCCGGG"_dna5};

    auto t = foo | std::views::take(1);                      // == [ [A,A,A,T,T,T] ]

    auto d = foo | seqan3::views::deep{std::views::take(1)}; // == [ [A], [C] ]
    // constructor arguments passed via {} and arguments to underlying view hardcoded inside

    // or with an alias defined previously
    auto e = foo | my::deep_take1;                           // == [ [A], [C] ]
}
