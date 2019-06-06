#include <vector>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/fm_index/all.hpp>

using namespace seqan3;

int main()
{

// TODO(comments): outputs [A, A, G] instead of AAG

//! [cycle]
std::vector<dna4> genome{"AATAATAAC"_dna4};
fm_index index{genome};                                 // build the index

auto cur = index.begin();                               // create a cursor
// cur.cycle_back();                                    // cycle_back on begin() is undefined behaviour!
cur.extend_right("AAC"_dna4);                           // search the sequence "AAC"
debug_stream << cur.path_label(genome) << '\n';         // outputs "AAC"
debug_stream << cur.last_rank() << '\n';                // outputs 1

cur.cycle_back();                                       // search the sequence "AAT"
debug_stream << cur.path_label(genome) << '\n';         // outputs "AAT"
debug_stream << cur.last_rank() << '\n';                // outputs 3

cur.cycle_back();                                       // "cur" doesn't change because the rightmost char is already the largest dna4 char.
debug_stream << cur.path_label(genome) << '\n';         // outputs "AAT"
debug_stream << cur.last_rank() << '\n';                // outputs 3
//! [cycle]

return 0;
}
