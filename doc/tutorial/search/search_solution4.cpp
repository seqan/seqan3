#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/algorithm/search.hpp>

using namespace seqan3;

int main()
{
    dna4_vector text{"CGCTGTCTGAAGGATGAGTGTCAGCCAGTGTAACCCGATGAGCTACCCAGTAGTCGAACTGGGCCAGACAACCCGGCGCTAATGCACTCA"_dna4};
    dna4_vector query{"GCT"_dna4};

    fm_index index{text};

    debug_stream << "Searching all hits\n";
    configuration const cfg_all = search_cfg::max_error{search_cfg::total{1}} |
                                  search_cfg::mode{search_cfg::all};
    auto results_all = search(query, index, cfg_all);
    debug_stream << "There are " << results_all.size() << " hits.\n";

    debug_stream << "Searching all best hits\n";
    configuration const cfg_all_best = search_cfg::max_error{search_cfg::total{1}} |
                                       search_cfg::mode{search_cfg::all_best};
    auto results_all_best = search(query, index, cfg_all_best);
    debug_stream << "There are " << results_all_best.size() << " hits.\n";

    debug_stream << "Searching best hit\n";
    configuration const cfg_best = search_cfg::max_error{search_cfg::total{1}} |
                                   search_cfg::mode{search_cfg::best};
    auto results_best = search(query, index, cfg_best);
    debug_stream << "There is " << results_best.size() << " hit.\n";

    debug_stream << "Searching all hits in the 1-stratum\n";
    configuration const cfg_strata = search_cfg::max_error{search_cfg::total{1}} |
                                     search_cfg::mode{search_cfg::strata{1}};
    auto results_strata = search(query, index, cfg_strata);
    debug_stream << "There are " << results_strata.size() << " hits.\n";
}
