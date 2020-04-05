#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/search.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>

using seqan3::operator""_dna4;

int main()
{
    seqan3::dna4_vector
                text{"CGCTGTCTGAAGGATGAGTGTCAGCCAGTGTAACCCGATGAGCTACCCAGTAGTCGAACTGGGCCAGACAACCCGGCGCTAATGCACTCA"_dna4};
    seqan3::dna4_vector query{"GCT"_dna4};

    seqan3::fm_index index{text};

    seqan3::debug_stream << "Searching all hits\n";
    seqan3::configuration const cfg_all = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1}} |
                                          seqan3::search_cfg::mode{seqan3::search_cfg::all};
    auto results_all = search(query, index, cfg_all);
    seqan3::debug_stream << "Hits: " << results_all << "\n";

    seqan3::debug_stream << "Searching all best hits\n";
    seqan3::configuration const cfg_all_best = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1}} |
                                               seqan3::search_cfg::mode{seqan3::search_cfg::all_best};
    auto results_all_best = search(query, index, cfg_all_best);
    seqan3::debug_stream << "Hits: " << results_all_best << "\n";

    seqan3::debug_stream << "Searching best hit\n";
    seqan3::configuration const cfg_best = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1}} |
                                           seqan3::search_cfg::mode{seqan3::search_cfg::best};
    auto results_best = search(query, index, cfg_best);
    seqan3::debug_stream << "Hits " << results_best << "\n";

    seqan3::debug_stream << "Searching all hits in the 1-stratum\n";
    seqan3::configuration const cfg_strata = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1}} |
                                             seqan3::search_cfg::mode{seqan3::search_cfg::strata{1}};
    auto results_strata = search(query, index, cfg_strata);
    seqan3::debug_stream << "Hits: " << results_strata << "\n";
}
