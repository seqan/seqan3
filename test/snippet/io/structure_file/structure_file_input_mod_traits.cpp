#include <sstream>

#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/structure_file/input.hpp>

// Define custom traits
struct my_traits : seqan3::structure_file_input_default_traits_rna
{
    using seq_alphabet = seqan3::rna4; // instead of rna5
};

auto input = R"(> S.cerevisiae_tRNA-PHE M10740/1-73
GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA
(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)
> example
UUGGAGUACACAACCUGUACACUCUUUC
..(((((..(((...)))..)))))... (-3.71))";

int main()
{
    seqan3::structure_file_input<my_traits> fin{std::istringstream{input}, seqan3::format_vienna{}};
}
