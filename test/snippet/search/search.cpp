#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/search.hpp>
#include <seqan3/search/fm_index/all.hpp>

int main()
{
    using seqan3::operator""_dna4;
    std::vector<seqan3::dna4_vector> genomes{"CGCTGTCTGAAGGATGAGTGTCAGCCAGTGTA"_dna4,
                                             "ACCCGATGAGCTACCCAGTAGTCGAACTG"_dna4,
                                             "GGCCAGACAACCCGGCGCTAATGCACTCA"_dna4};
    std::vector<seqan3::dna4_vector> queries{"GCT"_dna4, "ACCC"_dna4};

    // build an FM index
    seqan3::fm_index index{genomes};

    // search for the queries "GCT" and "ACCC"
    seqan3::debug_stream << seqan3::search(queries, index) << '\n';
    // This should result in: [[(0,1),(1,9),(2,16)],[(1,0),(1,12),(2,9)]], where the first list [(0,1),(1,9),(2,16)]
    // is the position list of "GCT" and the second one [(1,0),(1,12),(2,9)] the position list of "ACCC".
    // The first element of one tuple indicates in which sequence one of the queries is found
    // (eg. 0 stands for the first sequence in genomes). The second element of one tuple gives the position in a
    // sequence, where the query can be found. For example, (0,1) indicates that "GCT" can be found in the first
    // sequence at position 1, which is true for the here given example.

}
