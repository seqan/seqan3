// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2018, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::alignment_algorithm.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/view/drop_exactly.hpp>
#include <range/v3/view/zip.hpp>

#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/pairwise/policy/affine_gap_init_policy.hpp>
#include <seqan3/alignment/pairwise/policy/affine_gap_policy.hpp>
#include <seqan3/alignment/pairwise/policy/unbanded_dp_matrix_policy.hpp>
#include <seqan3/alignment/pairwise/select_result_type.hpp>

#include <seqan3/alignment/scoring/gap_scheme.hpp>
#include <seqan3/alignment/scoring/scoring_scheme_base.hpp>

#include <seqan3/core/metafunction/deferred_crtp_base.hpp>
#include <seqan3/io/stream/debug_stream.hpp>
#include <seqan3/range/view/take_exactly.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

/*!\brief A function object that computes pairwise alignments given two sequences.
 * \ingroup pairwise_alignment
 * \tparam config_t             The configuration type; must be of type seqan3::configuration.
 * \tparam algorithm_policies_t Template parameter pack with the policies to determine the execution of the algorithm;
 *                              must be wrapped as seqan3::detail::deferred_crtp_base.
 *
 * \details
 *
 * Computes the actual pairwise alignment given two sequences. The correct algorithm is configured while parsing
 * the seqan3::configuration within the detail::seqan3::configurator.
 */
template <typename config_t, typename ...algorithm_policies_t>
class alignment_algorithm :
    protected invoke_deferred_crtp_base<algorithm_policies_t, alignment_algorithm<algorithm_policies_t...>>...
{
public:
    /*!\name Constructor, destructor and assignment
     * \brief Defaulted all standard constructor.
     * \{
     */
    constexpr alignment_algorithm()                                        = default;
    constexpr alignment_algorithm(alignment_algorithm const &)             = default;
    constexpr alignment_algorithm(alignment_algorithm &&)                  = default;
    constexpr alignment_algorithm & operator=(alignment_algorithm const &) = default;
    constexpr alignment_algorithm & operator=(alignment_algorithm &&)      = default;
    ~alignment_algorithm()                                                 = default;

    /*!\brief Constructs the algorithm with the passed configuration.
     * \param cfg The configuration to be passed to the algorithm.
     *
     * \details
     *
     * The configuration is copied once to the heap during construction and maintained by a std::shared_ptr.
     * The configuration is not passed to the function-call-operator of this function object, in order to avoid
     * incompatible configurations between the passed configuration and the one used during configuration of this
     * class. Further will the function object be stored in a std::function which requires copyable objects and
     * in parallel executions the function object must be copied as well.
     */
    explicit constexpr alignment_algorithm(config_t const & cfg) : cfg_ptr{new config_t(cfg)}
    {}
    //!}

    /*!\brief Invokes the actual alignment computation given two sequences.
     * \tparam    first_batch_t  The type of the first sequence (or packed sequences); must model std::ForwardRange.
     * \tparam    second_batch_t The type of the second sequence (or packed sequences); must model std::ForwardRange.
     * \param[in] first_batch    The first sequence (or packed sequences).
     * \param[in] second_batch   The second sequence (or packed sequences).
     */
    template <std::ranges::ForwardRange first_batch_t, std::ranges::ForwardRange second_batch_t>
    auto operator()(first_batch_t const & first_batch, second_batch_t const & second_batch)
    {
        assert(cfg_ptr != nullptr);

        // We need to allocate the score_matrix and maybe the trace_matrix.
        this->allocate_score_matrix(first_batch, second_batch);

        // Initialize cache variables to keep frequently used variables close to the CPU registers.
        auto cache = this->setup_cache(get<align_cfg::gap>(*cfg_ptr).value);
        // auto cache = this->setup_cache(gap_scheme{gap_score{-1}, gap_open_score{-10}});

        initialize(cache);

        compute(first_batch, second_batch, cache);

        using result_t = typename select_result_type<first_batch_t, second_batch_t, config_t>::type;
        result_t res{};

        // Choose what needs to be computed.
        if constexpr (config_t::template exists<align_cfg::result<detail::with_score_type>>())
        {
            auto last_col = this->active_column();
            get<1>(res) = get<0>(*(std::ranges::end(last_col) - 1));
            return res;
        }
        return res;

    }

protected:

    /*!\brief Initializes the first column of the dynamic programming matrix.
     * \tparam         cache_t The cache type.
     * \param[in,out]  cache   The cache holding hot variables.
     */
    template <typename cache_t>
    void initialize(cache_t & cache)
    {
        auto col = this->active_column();

        this->init_origin_cell(*std::ranges::begin(col), cache);

        ranges::for_each(col | ranges::view::drop_exactly(1), [&cache, this](auto & cell)
        {
            this->init_column_cell(cell, cache);
        });
    }

    /*!\brief Compute the alignment by iterating over the dynamic programming matrix in a column wise manner.
     * \tparam        first_batch_t  The type of the first sequence (or packed sequences).
     * \tparam        second_batch_t The type of the second sequence (or packed sequences).
     * \tparam        cache_t        The type of the cache.
     * \param[in]     first_batch    The first sequence.
     * \param[in]     second_batch   The second sequence.
     * \param[in,out] cache          The cache holding hot variables.
     */
    template <typename first_batch_t,
              typename second_batch_t,
              typename cache_t>
    void compute(first_batch_t const & first_batch,
                 second_batch_t const & second_batch,
                 cache_t & cache)
    {
        auto const & score_scheme = get<align_cfg::scoring>(*cfg_ptr).value;
        ranges::for_each(first_batch, [&, this](auto seq1_value)
        {
            auto col = this->active_column();
            this->init_row_cell(*std::ranges::begin(col), cache);

            auto second_batch_it = std::ranges::begin(second_batch);
            ranges::for_each(col | ranges::view::drop_exactly(1), [&, this] (auto & cell)
            {
                this->compute_cell(cell, cache, score_scheme.score(seq1_value, *second_batch_it));
                ++second_batch_it;
            });
        });
    }

private:

    //!\brief The alignment configuration stored on the heap.
    std::shared_ptr<config_t> cfg_ptr{};
};

} // namespace seqan3::detail
