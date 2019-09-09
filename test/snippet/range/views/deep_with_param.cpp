#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/range/views/deep.hpp>
#include <seqan3/range/views/take.hpp>

namespace my
{
inline auto const deep_take = seqan3::views::deep{std::views::take};
}

int main()
{
    using seqan3::operator""_dna5;

    std::vector<seqan3::dna5_vector> foo{"AAATTT"_dna5, "CCCGGG"_dna5};

    auto t = foo | seqan3::views::take(1);                   // == [ [A,A,A,T,T,T] ]

    auto d = foo | seqan3::views::deep{std::views::take}(1); // == [ [A], [C] ]
    // constructor arguments passed via {} and arguments to underlying view passed via ()

    // In this case especially, an alias improves readability:
    auto e = foo | my::deep_take(1);                         // == [ [A], [C] ]
}
