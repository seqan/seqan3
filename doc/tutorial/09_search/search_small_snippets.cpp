// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include "cleanup.hpp"
seqan3::cleanup index_file{"index.file"};

#if SEQAN3_WITH_CEREAL
#    include <fstream>

#    include <cereal/archives/binary.hpp>
#endif //SEQAN3_WITH_CEREAL

#include <filesystem>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

using namespace std::string_literals;

int main()
{

    {
        //![text_collection]
        std::vector<std::string> text{{"Garfield the fat cat without a hat."},
                                      {"He is infinite, he is eternal."},
                                      {"Yet another text I have to think about."}};
        seqan3::fm_index index{text};
        seqan3::bi_fm_index bi_index{text};
        //![text_collection]
    }

#if SEQAN3_WITH_CEREAL
    {
//![store]
#    include <fstream> // for writing/reading files

#    include <cereal/archives/binary.hpp> // for storing/loading indices via cereal

        std::string text{"Garfield the fat cat without a hat."};
        seqan3::fm_index index{text};
        {
            std::ofstream os{"index.file", std::ios::binary};
            cereal::BinaryOutputArchive oarchive{os};
            oarchive(index);
        }
        //![store]
    }

    {
        //![load]
        // we need to tell the index that we work on a single text and a `char` alphabet before loading
        seqan3::fm_index<char, seqan3::text_layout::single> index;
        {
            std::ifstream is{"index.file", std::ios::binary};
            cereal::BinaryInputArchive iarchive{is};
            iarchive(index);
        }
        //![load]
    }
#endif //SEQAN3_WITH_CEREAL

    {
        //![error_search]
        std::string text{"Garfield the fat cat without a hat."};
        seqan3::fm_index index{text};
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{1}}
                                        | seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_count{0}}
                                        | seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_count{1}}
                                        | seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_count{1}};
        seqan3::debug_stream << search("cat"s, index, cfg) << '\n';
        // prints: [<query_id:0, reference_id:0, reference_pos:14>,
        //          <query_id:0, reference_id:0, reference_pos:17>,
        //          <query_id:0, reference_id:0, reference_pos:18>,
        //          <query_id:0, reference_id:0, reference_pos:32>]
        //![error_search]
    }

    {
        //![multiple_queries]
        std::string text{"Garfield the fat cat without a hat."};
        seqan3::fm_index index{text};
        std::vector<std::string> query{"cat"s, "hat"s};
        seqan3::debug_stream << search(query, index) << '\n';
        // prints: [<query_id:0, reference_id:0, reference_pos:17>,
        //          <query_id:1, reference_id:0, reference_pos:31>]
        //![multiple_queries]
    }

    {
        //![error_sum]
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{2}}
                                        | seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_count{2}}
                                        | seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_count{1}}
                                        | seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_count{1}};
        //![error_sum]
    }

    {
        //![hit_best]
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{1}}
                                        | seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_count{0}}
                                        | seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_count{1}}
                                        | seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_count{1}}
                                        | seqan3::search_cfg::hit_single_best{};
        //![hit_best]
    }

    {
        //![hit_strata]
        seqan3::configuration cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{2}}
                                  | seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_count{0}}
                                  | seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_count{1}}
                                  | seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_count{1}}
                                  | seqan3::search_cfg::hit_strata{2};
        using seqan3::get;                                    // Required to be able to find the correct get function.
        get<seqan3::search_cfg::hit_strata>(cfg).stratum = 1; // The stratum is now 1 and not 2 anymore.
        //![hit_strata]
    }

    //![hit_dynamic]
    seqan3::search_cfg::hit hit_dynamic{seqan3::search_cfg::hit_all{}}; // Initialise with hit_all configuration.

    bool hit_with_strata = static_cast<bool>(std::rand() & 1); // Either false or true.
    if (hit_with_strata)
        hit_dynamic = seqan3::search_cfg::hit_strata{2}; // Search instead with strata mode.

    seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{2}}
                                    | seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_count{0}}
                                    | seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_count{1}}
                                    | seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_count{1}}
                                    | hit_dynamic; // Build the configuration by adding the dynamic hit configuration.
    //![hit_dynamic]
}
