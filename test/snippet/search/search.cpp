#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/search.hpp>
#include <seqan3/search/fm_index/all.hpp>

int main()
{
    using namespace seqan3::literals;
    std::vector<seqan3::dna4_vector> genomes{"CGCTGTCTGAAGGATGAGTGTCAGCCAGTGTA"_dna4,
                                             "ACCCGATGAGCTACCCAGTAGTCGAACTG"_dna4,
                                             "GGCCAGACAACCCGGCGCTAATGCACTCA"_dna4};
    std::vector<seqan3::dna4_vector> queries{"GCT"_dna4, "ACCC"_dna4};

    // build an FM index
    seqan3::fm_index index{genomes};

    // search for the queries "GCT" and "ACCC"
    for (auto && result : seqan3::search(queries, index))
        seqan3::debug_stream << result << '\n';
    // This should result in:
    // <query_id:0, reference_id:0, reference_pos:1>
    // <query_id:0, reference_id:1, reference_pos:9>
    // <query_id:0, reference_id:2, reference_pos:16>
    // <query_id:1, reference_id:1, reference_pos:0>
    // <query_id:1, reference_id:1, reference_pos:12>
    // <query_id:1, reference_id:2, reference_pos:9>
}
