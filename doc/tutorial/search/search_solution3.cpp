#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/algorithm/search.hpp>
#include <seqan3/std/span>

using namespace seqan3;

int main()
{
    dna4_vector text{"CGCTGTCTGAAGGATGAGTGTCAGCCAGTGTAACCCGATGAGCTACCCAGTAGTCGAACTGGGCCAGACAACCCGGCGCTAATGCACTCA"_dna4};
    fm_index index{text};

    configuration const cfg = search_cfg::max_error{search_cfg::substitution{1}};

    auto results = search("GCT"_dna4, index, cfg);
    std::ranges::sort(results);

    debug_stream << "There are " << results.size() << " hits.\n";

    for (auto && pos : results)
    {
        debug_stream << "At position " << pos << ": "
                     << std::span{std::data(text) + pos, 3}
                     << '\n';
    }
}
