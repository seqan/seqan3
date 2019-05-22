#include <vector>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/fm_index/all.hpp>

using namespace seqan3;

int main()
{
    std::vector<std::vector<dna4>> genomes{"ATCTGACGAAGGCTAGCTAGCTAAGGGA"_dna4, "TAGCTGAAGCCATTGGCATCTGATCGGACT"_dna4,
                                           "ACTGAGCTCGTC"_dna4, "TGCATGCACCCATCGACTGACTG"_dna4, "GTACGTACGTTACG"_dna4};
    fm_index index{genomes};                                  // build the index

    auto it = index.begin();                                  // create an iterator
    it.extend_right("CTGA"_dna4);                             // search the pattern "CTGA"
    debug_stream << "Number of hits: " << it.count() << '\n'; // outputs: 5
    debug_stream << "Positions in the genomes: ";
    for (auto const & pos : it.locate())                      // outputs: (3,16) (2,1) (1,3) (0,2) (1,19)
        debug_stream << pos << ' ';
    debug_stream << '\n';
    return 0;
}
