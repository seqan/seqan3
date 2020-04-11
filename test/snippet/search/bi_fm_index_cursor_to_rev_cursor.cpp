#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>

int main()
{
    using seqan3::operator""_dna4;

    seqan3::debug_stream << "Example to_rev_cursor()\n";

    seqan3::dna4_vector genome{"GAATTAACGAAC"_dna4};
    seqan3::bi_fm_index index{genome};                              // build the bidirectional index

    auto cur = index.cursor();                                      // create a cursor
    cur.extend_right("AAC"_dna4);                                   // search the sequence "AAC"
    seqan3::debug_stream << cur.path_label(genome) << '\n';         // outputs "AAC"
    auto rev_cur = cur.to_rev_cursor();                             // unidirectional cursor on the text "CAAGCAATTAAG"
    auto genome_rev = genome | std::views::reverse;                 // create reverse text
    seqan3::debug_stream << rev_cur.path_label(genome_rev) << '\n'; // outputs "CAA"
    // Undefined behaviour! Cannot be called on the reversed cursor if the last extension on the bidirectional
    // cursor was to the right:
    // cur.cycle_back();
    // seqan3::debug_stream << cur.last_rank() << '\n';

    rev_cur.extend_right('G'_dna4);                                  // search the sequence "CAAG"
    seqan3::debug_stream << rev_cur.path_label(genome_rev) << '\n';  // outputs "CAAG"
    seqan3::debug_stream << rev_cur.last_rank() << '\n';             // outputs 2
    rev_cur.cycle_back();                                            // search the sequence "CAAT"
}
