#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>

int main()
{
    using seqan3::operator""_dna4;

    seqan3::debug_stream << "Example to_rev_cursor() on collections\n";

    std::vector<std::vector<seqan3::dna4>> genomes{"GAATTAACGAAC"_dna4, "TTTAACTTATC"_dna4};
    seqan3::bi_fm_index index{genomes};                  // build the bidirectional index

    auto cur = index.cursor();                           // create a cursor
    cur.extend_right("AAC"_dna4);                        // search the sequence "AAC"
    seqan3::debug_stream << cur.locate() << '\n';        // outputs [(0,9),(0,5),(1,3)]
    auto uni_it = cur.to_rev_cursor();                   // unidirectional cursor on the text "CTATTCAATTT|CAAGCAATTAAG"
    seqan3::debug_stream << uni_it.locate() << '\n';     // outputs [(1,4),(0,5),(1,0)] for "CAA"
}
