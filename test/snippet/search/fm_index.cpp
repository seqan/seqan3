#include <vector>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/fm_index/all.hpp>

using namespace seqan3;

int main()
{
std::vector<dna4> genome{"ATCGATCGAAGGCTAGCTAGCTAAGGGA"_dna4};
fm_index index{genome};                                    // build the index

auto cur = index.begin();                                  // create an cursor
cur.extend_right("AAGG"_dna4);                             // search the pattern "AAGG"
debug_stream << "Number of hits: " << cur.count() << '\n'; // outputs: 2
debug_stream << "Positions in the genome: ";
for (auto const & pos : cur.locate())                      // outputs: 8, 22
    debug_stream << pos << ' ';
debug_stream << '\n';
return 0;
}
