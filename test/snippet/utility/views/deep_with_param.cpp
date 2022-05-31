#include <ranges>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/views/deep.hpp>

namespace my
{
inline auto const deep_take = seqan3::views::deep{std::views::take};
}

int main()
{
    using namespace seqan3::literals;

    std::vector<seqan3::dna5_vector> sequences{"AAATTT"_dna5, "CCCGGG"_dna5};

    seqan3::debug_stream << (sequences | std::views::take(1)) << '\n'; // [A,A,A,T,T,T]

    seqan3::debug_stream << (sequences | seqan3::views::deep{std::views::take(1)}) << '\n'; // [A,C]
    // constructor arguments passed via {} and arguments to underlying view passed via ()

    // In this case especially, an alias improves readability:
    seqan3::debug_stream << (sequences | my::deep_take(1)) << '\n'; // [A,C]
}
