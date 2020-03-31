#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/search.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/std/span>

using seqan3::operator""_dna4;

void run_text_single()
{
    seqan3::dna4_vector
                text{"CGCTGTCTGAAGGATGAGTGTCAGCCAGTGTAACCCGATGAGCTACCCAGTAGTCGAACTGGGCCAGACAACCCGGCGCTAATGCACTCA"_dna4};
    seqan3::dna4_vector query{"GCT"_dna4};
    seqan3::fm_index index{text};

    seqan3::debug_stream << "Searching all best hits allowing for 1 error in a single text\n";

    seqan3::configuration const search_config = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1}} |
                                                seqan3::search_cfg::mode{seqan3::search_cfg::all_best};
    seqan3::configuration const align_config = seqan3::align_cfg::edit |
                                               seqan3::align_cfg::aligned_ends{seqan3::free_ends_first} |
                                               seqan3::align_cfg::result{seqan3::with_alignment};

    auto results = search(query, index, search_config);

    seqan3::debug_stream << "-----------------\n";

    for (auto && pos : results)
    {
        size_t start = pos.second ? pos.second - 1 : 0;
        std::span text_view{std::data(text) + start, query.size() + 1};

        for (auto && res : align_pairwise(std::tie(text_view, query), align_config))
        {
            auto && [aligned_database, aligned_query] = res.alignment();
            seqan3::debug_stream << "score:    " << res.score() << '\n';
            seqan3::debug_stream << "database: " << aligned_database << '\n';
            seqan3::debug_stream << "query:    "  << aligned_query << '\n';
            seqan3::debug_stream << "=============\n";
        }
    }
}

void run_text_collection()
{
    std::vector<seqan3::dna4_vector> text{"CGCTGTCTGAAGGATGAGTGTCAGCCAGTGTA"_dna4,
                                          "ACCCGATGAGCTACCCAGTAGTCGAACTG"_dna4,
                                          "GGCCAGACAACCCGGCGCTAATGCACTCA"_dna4};
    seqan3::dna4_vector query{"GCT"_dna4};
    seqan3::fm_index index{text};

    seqan3::debug_stream << "Searching all best hits allowing for 1 error in a text collection\n";

    seqan3::configuration const search_config = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1}} |
                                                seqan3::search_cfg::mode{seqan3::search_cfg::all_best};
    seqan3::configuration const align_config = seqan3::align_cfg::edit |
                                               seqan3::align_cfg::aligned_ends{seqan3::free_ends_first} |
                                               seqan3::align_cfg::result{seqan3::with_alignment};

    seqan3::debug_stream << "-----------------\n";

    for (auto [query_idx, text_pos] : search(query, index, search_config))
    {
        size_t start = text_pos.second ? text_pos.second - 1 : 0;
        std::span text_view{std::data(text[text_pos.first]) + start, query.size() + 1};

        for (auto && res : align_pairwise(std::tie(text_view, query), align_config))
        {
            auto && [aligned_database, aligned_query] = res.alignment();
            seqan3::debug_stream << "score:    " << res.score() << '\n';
            seqan3::debug_stream << "database: " << aligned_database << '\n';
            seqan3::debug_stream << "query:    "  << aligned_query << '\n';
            seqan3::debug_stream << "=============\n";
        }
    }
}

int main()
{
   run_text_single();
   seqan3::debug_stream << '\n';
   run_text_collection();
}
