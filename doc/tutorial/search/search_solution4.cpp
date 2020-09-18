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
    seqan3::configuration const cfg_all = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{1}} |
                                          seqan3::search_cfg::hit_all{};
    auto results_all = search(query, index, cfg_all);
    // Attention: results_all is a pure std::ranges::input_range,
    //            so after calling std::ranges::distance, you cannot iterate over it again!
    seqan3::debug_stream << "There are " << std::ranges::distance(results_all) << " hits.\n";

    seqan3::debug_stream << "Searching all best hits\n";
    seqan3::configuration const cfg_all_best = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{1}} |
                                               seqan3::search_cfg::hit_all_best{};
    auto results_all_best = search(query, index, cfg_all_best);
    // Attention: results_all_best is a pure std::ranges::input_range,
    //            so after calling std::ranges::distance, you cannot iterate over it again!
    seqan3::debug_stream << "There are " << std::ranges::distance(results_all_best) << " hits.\n";

    seqan3::debug_stream << "Searching best hit\n";
    seqan3::configuration const cfg_best = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{1}} |
                                           seqan3::search_cfg::hit_single_best;
    auto results_best = search(query, index, cfg_best);
    // Attention: results_best is a pure std::ranges::input_range,
    //            so after calling std::ranges::distance, you cannot iterate over it again!
    seqan3::debug_stream << "There is " << std::ranges::distance(results_best) << " hits.\n";

    seqan3::debug_stream << "Searching all hits in the 1-stratum\n";
    seqan3::configuration const cfg_strata = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{1}} |
                                             seqan3::search_cfg::hit_strata{1};
    auto results_strata = search(query, index, cfg_strata);
    // Attention: results_strata is a pure std::ranges::input_range,
    //            so after calling std::ranges::distance, you cannot iterate over it again!
    seqan3::debug_stream << "There are " << std::ranges::distance(results_strata) << " hits.\n";
}
