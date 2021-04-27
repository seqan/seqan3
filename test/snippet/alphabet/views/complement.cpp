#include <seqan3/std/ranges>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/complement.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dna5_vector foo{"ACGTA"_dna5};

    // pipe notation
    auto v = foo | seqan3::views::complement;                        // == "TGCAT"

    // function notation
    auto v2(seqan3::views::complement(foo));                         // == "TGCAT"

    // generate the reverse complement:
    auto v3 = foo | seqan3::views::complement | std::views::reverse; // == "TACGT"
}
