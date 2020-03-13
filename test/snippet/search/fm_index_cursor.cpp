#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/fm_index/all.hpp>

int main()
{
    using seqan3::operator""_dna4;

    std::vector<seqan3::dna4> genome{"AATAATAAC"_dna4};
    seqan3::fm_index index{genome};                                 // build the index

    auto cur = index.cursor();                                      // create a cursor
    // cur.cycle_back();                                            // cycle_back on begin() is undefined behaviour!
    cur.extend_right("AAC"_dna4);                                   // search the sequence "AAC"
    seqan3::debug_stream << cur.path_label(genome) << '\n';         // outputs "AAC"
    seqan3::debug_stream << cur.last_rank() << '\n';                // outputs 1

    cur.cycle_back();                                               // search the sequence "AAT"
    seqan3::debug_stream << cur.path_label(genome) << '\n';         // outputs "AAT"
    seqan3::debug_stream << cur.last_rank() << '\n';                // outputs 3

    cur.cycle_back();                                               // "cur" doesn't change because the rightmost char
                                                                    // is already the largest dna4 char.
    seqan3::debug_stream << cur.path_label(genome) << '\n';         // outputs "AAT"
    seqan3::debug_stream << cur.last_rank() << '\n';                // outputs 3
}
