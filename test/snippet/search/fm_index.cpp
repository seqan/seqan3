#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/fm_index/all.hpp>

int main()
{
    using namespace seqan3::literals;

    std::vector<seqan3::dna4> genome{"ATCGATCGAAGGCTAGCTAGCTAAGGGA"_dna4};
    seqan3::fm_index index{genome};                                    // build the index

    auto cur = index.cursor();                                         // create a cursor
    cur.extend_right("AAGG"_dna4);                                     // search the pattern "AAGG"
    seqan3::debug_stream << "Number of hits: " << cur.count() << '\n'; // outputs: 2
    seqan3::debug_stream << "Positions in the genome: ";
    for (auto && pos : cur.locate())                                   // outputs: (0, 8), (0, 22)
        seqan3::debug_stream << pos << ' ';
    seqan3::debug_stream << '\n';
    return 0;
}
