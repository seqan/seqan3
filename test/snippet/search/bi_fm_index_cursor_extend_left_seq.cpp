#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>

int main()
{
    using seqan3::operator""_dna4;

    seqan3::debug_stream << "Example extend_left(seq)\n";

    seqan3::dna4_vector genome{"GAATTAATGAAC"_dna4};
    seqan3::bi_fm_index index{genome};                      // build the bidirectional index

    auto cur = index.cursor();                              // create a cursor
    cur.extend_right("AAC"_dna4);                           // search the sequence "AAC"
    seqan3::debug_stream << cur.path_label(genome) << '\n'; // outputs "AAC"
    cur.extend_left("ATG"_dna4);                            // extend the query to "ATGAAC"
                                                            // The rightmost character of "ATG" is extended to the left first.
    seqan3::debug_stream << cur.path_label(genome) << '\n'; // outputs "ATGAAC"
}
