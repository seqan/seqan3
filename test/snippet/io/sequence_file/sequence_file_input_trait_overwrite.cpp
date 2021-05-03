#include <sstream>

#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sequence_file/input.hpp>

auto input = R"(> TEST1
ACGT
> Test2
AGGCTGA
> Test3
GGAGTATAATATATATATATATAT)";

struct my_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4; // instead of dna5

    template <typename alph>
    using sequence_container = seqan3::bitpacked_sequence<alph>; // must be defined as a template!
};

int main()
{
    seqan3::sequence_file_input<my_traits> fin{std::istringstream{input}, seqan3::format_fasta{}};
}
