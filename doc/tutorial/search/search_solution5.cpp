#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/algorithm/search.hpp>
#include <seqan3/std/span>

using namespace seqan3;

void run_text_single()
{
    dna4_vector text{"CGCTGTCTGAAGGATGAGTGTCAGCCAGTGTAACCCGATGAGCTACCCAGTAGTCGAACTGGGCCAGACAACCCGGCGCTAATGCACTCA"_dna4};
    dna4_vector query{"GCT"_dna4};
    fm_index index{text};

    debug_stream << "Searching all best hits allowing for 1 error in a single text\n";

    configuration const search_config = search_cfg::max_error{search_cfg::total{1}} |
                                        search_cfg::mode{search_cfg::all_best};
    configuration const align_config = align_cfg::edit |
                                       align_cfg::aligned_ends{free_ends_first} |
                                       align_cfg::result{with_alignment};

    auto results = search(query, index, search_config);

    debug_stream << "There are " << results.size() << " hits.\n";
    debug_stream << "-----------------\n";

    for (auto && pos : results)
    {
        size_t start = pos ? pos - 1 : 0;
        std::span text_view{std::data(text) + start, query.size() + 1};

        for (auto && res : align_pairwise(std::tie(text_view, query), align_config))
        {
            auto && [aligned_database, aligned_query] = res.alignment();
            debug_stream << "score:    " << res.score() << '\n';
            debug_stream << "database: " << aligned_database << '\n';
            debug_stream << "query:    "  << aligned_query << '\n';
            debug_stream << "=============\n";
        }
    }
}

void run_text_collection()
{
    std::vector<dna4_vector> text{"CGCTGTCTGAAGGATGAGTGTCAGCCAGTGTA"_dna4,
                                  "ACCCGATGAGCTACCCAGTAGTCGAACTG"_dna4,
                                  "GGCCAGACAACCCGGCGCTAATGCACTCA"_dna4};
    dna4_vector query{"GCT"_dna4};
    fm_index index{text};

    debug_stream << "Searching all best hits allowing for 1 error in a text collection\n";

    configuration const search_config = search_cfg::max_error{search_cfg::total{1}} |
                                        search_cfg::mode{search_cfg::all_best};
    configuration const align_config = align_cfg::edit |
                                       align_cfg::aligned_ends{free_ends_first} |
                                       align_cfg::result{with_alignment};

    auto results = search(query, index, search_config);

    debug_stream << "There are " << results.size() << " hits.\n";
    debug_stream << "-----------------\n";

    for (auto [idx, pos] : results)
    {
        size_t start = pos ? pos - 1 : 0;
        std::span text_view{std::data(text[idx]) + start, query.size() + 1};

        for (auto && res : align_pairwise(std::tie(text_view, query), align_config))
        {
            auto && [aligned_database, aligned_query] = res.alignment();
            debug_stream << "score:    " << res.score() << '\n';
            debug_stream << "database: " << aligned_database << '\n';
            debug_stream << "query:    "  << aligned_query << '\n';
            debug_stream << "=============\n";
        }
    }
}

int main()
{
   run_text_single();
   debug_stream << '\n';
   run_text_collection();
}
