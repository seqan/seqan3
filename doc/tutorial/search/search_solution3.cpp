#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/search.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/std/span>

using seqan3::operator""_dna4;

int main()
{
    seqan3::dna4_vector
                text{"CGCTGTCTGAAGGATGAGTGTCAGCCAGTGTAACCCGATGAGCTACCCAGTAGTCGAACTGGGCCAGACAACCCGGCGCTAATGCACTCA"_dna4};
    seqan3::fm_index index{text};

    seqan3::configuration const cfg = seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_count{1}};

    for (auto && result : search("GCT"_dna4, index, cfg))
    {
        seqan3::debug_stream << "At position " << result.reference_begin_position() << ": "
                             << std::span{std::data(text) + result.reference_begin_position(), 3}
                             << '\n';
    }
}
