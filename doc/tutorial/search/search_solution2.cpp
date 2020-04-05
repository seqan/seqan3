#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/search.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>

using seqan3::operator""_dna4;

void run_text_single()
{
    seqan3::dna4_vector
                text{"CGCTGTCTGAAGGATGAGTGTCAGCCAGTGTAACCCGATGAGCTACCCAGTAGTCGAACTGGGCCAGACAACCCGGCGCTAATGCACTCA"_dna4};
    seqan3::fm_index index{text};

    seqan3::debug_stream << "=====   Running on a single text   =====\n";
    seqan3::debug_stream << "The positions are " << search("GCT"_dna4, index) << '\n';
}

void run_text_collection()
{
    std::vector<seqan3::dna4_vector> text{"CGCTGTCTGAAGGATGAGTGTCAGCCAGTGTA"_dna4,
                                          "ACCCGATGAGCTACCCAGTAGTCGAACTG"_dna4,
                                          "GGCCAGACAACCCGGCGCTAATGCACTCA"_dna4};
    seqan3::fm_index index{text};

    seqan3::debug_stream << "===== Running on a text collection =====\n";
    seqan3::debug_stream << "The positions are " << search("GCT"_dna4, index) << '\n';
}

int main()
{
   run_text_single();
   seqan3::debug_stream << '\n';
   run_text_collection();
}
