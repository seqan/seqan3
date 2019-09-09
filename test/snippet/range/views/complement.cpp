#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/range/views/complement.hpp>
#include <seqan3/std/ranges>

int main()
{
    using seqan3::operator""_dna5;

    seqan3::dna5_vector foo{"ACGTA"_dna5};

    // pipe notation
    auto v = foo | seqan3::views::complement;                        // == "TGCAT"

    // function notation
    auto v2(seqan3::views::complement(foo));                         // == "TGCAT"

    // generate the reverse complement:
    auto v3 = foo | seqan3::views::complement | std::views::reverse; // == "TACGT"
}
