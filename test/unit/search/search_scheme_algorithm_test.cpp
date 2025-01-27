// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <algorithm>
#include <type_traits>

#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream/range.hpp>
#include <seqan3/search/configuration/default_configuration.hpp>
#include <seqan3/search/detail/policy_max_error.hpp>
#include <seqan3/search/detail/policy_search_result_builder.hpp>
#include <seqan3/search/detail/search_configurator.hpp>
#include <seqan3/search/detail/search_scheme_algorithm.hpp>
#include <seqan3/search/detail/unidirectional_search_algorithm.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/utility/range/to.hpp>
#include <seqan3/utility/views/slice.hpp>

#include "helper.hpp"
#include "helper_search_scheme.hpp"

// Uses the trivial search of the unidirectional search algorithm.
// The algorithm is configured with the corrsponding configuration types.
// To modify the trivial search use the configuration settings of the search algorithm.
template <typename index_t, typename query_t, typename delegate_t>
static void search_trivial(index_t const & index,
                           query_t const & query,
                           seqan3::detail::search_param error_left,
                           delegate_t && delegate)
{
    using namespace seqan3::detail;

    // Configure the algorithm according to the given specifications.
    auto cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{error_left.total}}
             | seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_count{error_left.substitution}}
             | seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_count{error_left.insertion}}
             | seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_count{error_left.deletion}}
             | seqan3::search_cfg::hit_all{} | seqan3::search_cfg::output_index_cursor{};

    auto indexed_query = std::pair{size_t{0}, query};
    auto algo =
        std::get<0>(seqan3::detail::search_configurator::configure_algorithm<decltype(indexed_query)>(cfg, index));

    // Call the algorithm and call the delegate with the returned index cursor.
    algo(indexed_query,
         [&](auto result)
         {
             delegate(result.index_cursor());
         });
}

template <typename text_t>
inline void test_search_hamming(auto index,
                                text_t const & text,
                                auto const & search,
                                uint64_t const query_length,
                                std::vector<uint8_t> const & error_distribution,
                                size_t const seed,
                                auto const & blocks_length,
                                auto const & ordered_blocks_length,
                                uint64_t const start_pos)
{
    using char_t = typename text_t::value_type;

    uint64_t const pos = std::rand() % (text.size() - query_length + 1);
    text_t const orig_query = text | seqan3::views::slice(pos, pos + query_length) | seqan3::ranges::to<text_t>();

    // Modify query s.t. it has errors matching error_distribution.
    auto query = orig_query;
    auto it = index.cursor();
    uint64_t current_blocks_length = 0;
    for (uint8_t block = 0; block < search.blocks(); ++block)
    {
        uint64_t const single_block_length = ordered_blocks_length[block];
        EXPECT_LE(error_distribution[block], single_block_length);
        if (error_distribution[block] > single_block_length)
        {
            seqan3::debug_stream_type<char> cerr{std::cerr};
            cerr << "Error in block " << block << "(+ 1): " << error_distribution[block]
                 << " errors cannot fit into a block of length " << single_block_length << "." << '\n'
                 << "Error Distribution: " << error_distribution << '\n';
            exit(1);
        }

        // Choose random positions in the query sequence for substitutions. Repeat until all error positions are unique.
        std::vector<uint8_t> error_positions(error_distribution[block]);
        do
        {
            error_positions.clear();
            for (uint8_t error = 0; error < error_distribution[block]; ++error)
                error_positions.push_back(std::rand() % single_block_length);
            std::sort(error_positions.begin(), error_positions.end());
        }
        while (std::adjacent_find(error_positions.begin(), error_positions.end()) != error_positions.end());

        // Construct query sequence with chosen error positions.
        for (size_t error = 0; error < error_positions.size(); ++error)
        {
            uint64_t const pos = error_positions[error] + current_blocks_length;
            // Decrease alphabet size by one because we don't want to replace query[pos], with the same character.
            uint8_t new_rank = std::rand() % (seqan3::alphabet_size<char_t> - 1);
            // If it is a match now, it can't be the highest rank of the alphabet. Thus we set it to the highest rank.
            if (new_rank == seqan3::to_rank(query[pos]))
                new_rank = seqan3::alphabet_size<char_t> - 1;
            seqan3::assign_rank_to(new_rank, query[pos]);
        }
        current_blocks_length += single_block_length;
    }

    std::vector<uint64_t> hits_trivial, hits_ss;

    auto delegate_trivial = [&hits_trivial](auto const & it)
    {
        for (auto && res : it.locate())
            hits_trivial.push_back(res.second);
    };

    auto delegate_ss = [&hits_ss](auto const & it)
    {
        for (auto && res : it.locate())
            hits_ss.push_back(res.second);
    };

    auto remove_predicate_ss = [&text, &orig_query, query_length](uint64_t const hit)
    {
        seqan3::dna4_vector matched_seq =
            text | seqan3::views::slice(hit, hit + query_length) | seqan3::ranges::to<seqan3::dna4_vector>();
        return (matched_seq != orig_query);
    };

    auto remove_predicate_trivial = [&](uint64_t const hit)
    {
        // filter only correct error distributions
        seqan3::dna4_vector matched_seq =
            text | seqan3::views::slice(hit, hit + query_length) | seqan3::ranges::to<seqan3::dna4_vector>();

        if (orig_query != matched_seq)
            return true;

        uint64_t lb = 0, rb = 0;
        for (uint8_t block = 0; block < search.blocks(); ++block)
        {
            rb += ordered_blocks_length[block];

            uint8_t errors = 0;
            for (uint64_t i = lb; i < rb; ++i)
                if (hit + i >= text.size())
                    ++errors;
                else
                    errors += query[i] != text[hit + i];
            if (errors != error_distribution[block])
                return true;
            lb += ordered_blocks_length[block];
        }
        return false;
    };

    uint8_t const total = search.u.back();
    uint8_t const substitution = std::rand() % (total + 1);

    seqan3::detail::search_param error_left{total, substitution, 0, 0};

    // Find all hits using search schemes.
    seqan3::detail::search_ss<
        false>(it, query, start_pos, start_pos + 1, 0, 0, true, search, blocks_length, error_left, delegate_ss);

    // Find all hits using trivial backtracking.
    search_trivial(index, query, error_left, delegate_trivial);

    // Eliminate hits that we are not interested in (based on the search and chosen error distribution)
    hits_ss.erase(std::remove_if(hits_ss.begin(), hits_ss.end(), remove_predicate_ss), hits_ss.end());

    hits_trivial.erase(std::remove_if(hits_trivial.begin(), hits_trivial.end(), remove_predicate_trivial),
                       hits_trivial.end());

    // Eliminate duplicates
    hits_ss = seqan3::uniquify(hits_ss);
    hits_trivial = seqan3::uniquify(hits_trivial);

    EXPECT_EQ(hits_ss, hits_trivial);
    if (hits_ss != hits_trivial)
    {
        seqan3::debug_stream_type<char> cerr{std::cerr};
        cerr << "Seed: " << seed << '\n'
             << "Text: " << text << '\n'
             << "Query: " << query << '\n'
             << "Errors: " << total << ", " << substitution << '\n'
             << "SS hits: " << hits_ss << '\n'
             << "Trivial hits: " << hits_trivial << '\n';
    }
}

template <typename search_scheme_t>
inline void
test_search_scheme_hamming(search_scheme_t const & search_scheme, size_t const seed, uint64_t const iterations)
{
    seqan3::dna4_vector text;

    search_scheme_t ordered_search_scheme;
    std::vector<std::vector<std::vector<uint8_t>>> error_distributions(search_scheme.size());

    // Calculate all error distributions and sort each of them (from left to right).
    uint8_t max_error = 0;
    for (size_t search_id = 0; search_id < search_scheme.size(); ++search_id)
    {
        ordered_search_scheme[search_id] = search_scheme[search_id];
        seqan3::search_error_distribution(error_distributions[search_id], search_scheme[search_id]);
        for (std::vector<uint8_t> & resElem : error_distributions[search_id])
            seqan3::order_search_vector(resElem, search_scheme[search_id]);
        max_error = std::max(max_error, search_scheme[search_id].u.back());
    }

    for (uint64_t text_length = 10; text_length < 10000; text_length *= 10)
    {
        uint64_t const query_length_min = std::max<uint64_t>(3, search_scheme.front().blocks() * max_error);
        uint64_t const query_length_max = std::min<uint64_t>(16, text_length);

        text = seqan3::test::generate_sequence<seqan3::dna4>(text_length, 0 /*variance*/, seed);
        seqan3::bi_fm_index index(text);

        for (uint64_t i = 0; i < iterations; ++i)
        {
            for (uint64_t query_length = query_length_min; query_length < query_length_max; ++query_length)
            {
                auto const block_info = search_scheme_block_info(search_scheme, query_length);
                for (size_t search_id = 0; search_id < search_scheme.size(); ++search_id)
                {
                    auto const & [blocks_length, start_pos] = block_info[search_id];

                    std::vector<uint64_t> ordered_blocks_length;
                    seqan3::get_ordered_search(search_scheme[search_id],
                                               blocks_length,
                                               ordered_search_scheme[search_id],
                                               ordered_blocks_length);

                    for (auto && error_distribution : error_distributions[search_id])
                    {
                        test_search_hamming(index,
                                            text,
                                            search_scheme[search_id],
                                            query_length,
                                            error_distribution,
                                            seed,
                                            blocks_length,
                                            ordered_blocks_length,
                                            start_pos);
                    }
                }
            }
        }
    }
}

template <typename search_scheme_t>
inline void test_search_scheme_edit(search_scheme_t const & search_scheme, size_t const seed, uint64_t const iterations)
{
    seqan3::dna4_vector text, query;

    // retrieve maximum number of errors from search_scheme
    uint8_t max_error = 0;
    for (auto const & search : search_scheme)
        max_error = std::max(max_error, search.u.back());

    for (uint64_t text_length = 10; text_length < 10000; text_length *= 10)
    {
        uint64_t const query_length_min = std::max<uint64_t>(3, search_scheme.front().blocks() * max_error);
        uint64_t const query_length_max = std::min<uint64_t>(16, text_length);

        text = seqan3::test::generate_sequence<seqan3::dna4>(text_length, 0 /*variance*/, seed);
        seqan3::bi_fm_index index(text);

        uint8_t const substitution = std::rand() % (max_error + 1);
        uint8_t const insertion = std::rand() % (max_error + 1);
        uint8_t const deletion = std::rand() % (max_error + 1);
        seqan3::detail::search_param error_left{max_error, substitution, insertion, deletion};

        for (uint64_t i = 0; i < iterations; ++i)
        {
            for (uint64_t query_length = query_length_min; query_length < query_length_max; ++query_length)
            {
                query = seqan3::test::generate_sequence<seqan3::dna4>(query_length, 0 /*variance*/, seed);

                std::vector<uint64_t> hits_trivial, hits_ss;

                auto delegate_trivial = [&hits_trivial](auto const & it)
                {
                    for (auto && res : it.locate())
                        hits_trivial.push_back(res.second);
                };

                auto delegate_ss = [&hits_ss](auto const & it)
                {
                    for (auto && res : it.locate())
                        hits_ss.push_back(res.second);
                };

                // Find all hits using search schemes.
                seqan3::detail::search_ss<false>(index, query, error_left, search_scheme, delegate_ss);
                // Find all hits using trivial backtracking.
                search_trivial(index, query, error_left, delegate_trivial);

                // Eliminate duplicates
                hits_ss = seqan3::uniquify(hits_ss);
                hits_trivial = seqan3::uniquify(hits_trivial);

                EXPECT_EQ(hits_ss, hits_trivial);
                if (hits_ss != hits_trivial)
                {
                    seqan3::debug_stream_type<char> cerr{std::cerr};
                    cerr << "Seed: " << seed << '\n'
                         << "Text: " << text << '\n'
                         << "Query: " << query << '\n'
                         << "Errors: " << max_error << ", " << substitution << ", " << insertion << ", " << deletion
                         << '\n'
                         << "SS hits: " << hits_ss << '\n'
                         << "Trivial hits: " << hits_trivial << '\n';
                }
            }
        }
    }
}

TEST(search_scheme_test, search_scheme_hamming)
{
    size_t seed = 42;

    test_search_scheme_hamming(seqan3::detail::optimum_search_scheme<0, 0>, seed, 10);
    test_search_scheme_hamming(seqan3::detail::optimum_search_scheme<0, 1>, seed, 10);
    test_search_scheme_hamming(seqan3::detail::optimum_search_scheme<1, 1>, seed, 10);
    test_search_scheme_hamming(seqan3::detail::optimum_search_scheme<0, 2>, seed, 10);
    test_search_scheme_hamming(seqan3::detail::optimum_search_scheme<1, 2>, seed, 10);
    test_search_scheme_hamming(seqan3::detail::optimum_search_scheme<2, 2>, seed, 10);
    test_search_scheme_hamming(seqan3::detail::optimum_search_scheme<0, 3>, seed, 10);
    test_search_scheme_hamming(seqan3::detail::optimum_search_scheme<1, 3>, seed, 10);
    test_search_scheme_hamming(seqan3::detail::optimum_search_scheme<2, 3>, seed, 10);
    // test_search_scheme_hamming(seqan3::detail::optimum_search_scheme<3, 3>, seed, 10);
}

TEST(search_scheme_test, search_scheme_edit)
{
    size_t seed = 42;

    // TODO: test with lower bounds != 0.
    // For that we need alignment statistics to know the number of errors spent in search_trivial
    test_search_scheme_edit(seqan3::detail::optimum_search_scheme<0, 0>, seed, 10);
    test_search_scheme_edit(seqan3::detail::optimum_search_scheme<0, 1>, seed, 10);
    test_search_scheme_edit(seqan3::detail::optimum_search_scheme<0, 2>, seed, 10);
    test_search_scheme_edit(seqan3::detail::optimum_search_scheme<0, 3>, seed, 10);
}
