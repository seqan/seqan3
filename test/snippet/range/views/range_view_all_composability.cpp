#include <ranges>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dna4_vector vec{"ACGGTC"_dna4};
    // views can be composed iteratively
    auto vec_view3 = vec | std::views::reverse;
    auto vec_view4 = vec_view3 | seqan3::views::complement;

    // or in one line similar to the unix command line
    auto vec_view5 = vec | seqan3::views::complement | std::views::reverse;

    // vec_view4 and vec_view5 are the reverse complement of "ACGGTC": "GACCGT"
    seqan3::debug_stream << vec_view4 << '\n'; // GACCGT
    seqan3::debug_stream << vec_view5 << '\n'; // GACCGT
}
