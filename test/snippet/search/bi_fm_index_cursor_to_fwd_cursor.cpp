#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>

int main()
{
    using seqan3::operator""_dna4;

    seqan3::debug_stream << "Example to_fwd_cursor()\n";

    seqan3::dna4_vector genome{"GAATTAACGAAC"_dna4};
    seqan3::bi_fm_index index{genome};                          // build the bidirectional index

    auto cur = index.cursor();                                  // create a cursor
    cur.extend_left("AAC"_dna4);                                // search the sequence "AAC"
    seqan3::debug_stream << cur.path_label(genome) << '\n';     // outputs "AAC"
    auto uni_it = cur.to_fwd_cursor();                          // unidirectional cursor on the text "GAATTAACGAAC"
    seqan3::debug_stream << uni_it.path_label(genome) << '\n';  // outputs "AAC"
    // Undefined behaviour! Cannot be called on the forward cursor if the last extension on the bidirectional
    // cursor was to the left:
    // cur.cycle_back();
    // seqan3::debug_stream << cur.last_rank() << '\n';

    uni_it.extend_right('G'_dna4);                              // search the sequence "AACG"
    seqan3::debug_stream << uni_it.path_label(genome) << '\n';  // outputs "AACG"
    seqan3::debug_stream << uni_it.last_rank() << '\n';         // outputs 2
    uni_it.cycle_back();                                        // returns false since there is no sequence "AACT" in the text.
}
