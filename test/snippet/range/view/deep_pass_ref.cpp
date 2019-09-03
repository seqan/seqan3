#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/range/view/deep.hpp>
#include <seqan3/range/view/take.hpp>

namespace my
{
inline auto const deep_take = seqan3::view::deep{std::view::take};
}

int main()
{
    using seqan3::operator""_dna5;

    std::vector<seqan3::dna5_vector> foo{"AAATTT"_dna5, "CCCGGG"_dna5};

    int i = 3;
    auto f = foo | my::deep_take(i); // takes `i` as a reference!
}
