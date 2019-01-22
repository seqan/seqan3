// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::alignment_selector.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <functional>
#include <tuple>
#include <utility>
#include <vector>

#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/pairwise/policy/all.hpp>
#include <seqan3/alignment/pairwise/alignment_algorithm.hpp>
#include <seqan3/alignment/pairwise/align_result_selector.hpp>
#include <seqan3/alignment/pairwise/align_result.hpp>
#include <seqan3/alignment/pairwise/edit_distance_unbanded.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/core/metafunction/deferred_crtp_base.hpp>
#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/core/metafunction/template_inspection.hpp>
#include <seqan3/core/type_list.hpp>
#include <seqan3/range/view/persist.hpp>

namespace seqan3::detail
{

/*!\brief Provides several contracts to test when configuring the alignment algorithm.
 * \ingroup pairwise_alignment
 * \tparam range_type            The type of the range containing sequences to be aligned.
 * \tparam alignment_config_type The type of the alignment configuration.
 *
 * \details
 *
 * This stateless helper class provides several contract testing functions for the alignment configuration.
 */
template <typename range_type,
          typename alignment_config_type>
struct alignment_contract
{
private:
    /*!\brief Auxiliary member types
     * \{
     */
    using first_seq_t  = std::tuple_element_t<0, value_type_t<std::ranges::iterator_t<range_type>>>;
    using second_seq_t = std::tuple_element_t<1, value_type_t<std::ranges::iterator_t<range_type>>>;
    //!\}

public:
    //!\brief Tests whether the value type of `range_type` is a tuple with exactly 2 members.
    constexpr static bool expects_tuple_like_value_type()
    {
        return tuple_like_concept<alignment_config_type> &&
               std::tuple_size_v<value_type_t<std::ranges::iterator_t<range_type>>> == 2;
    }

    //!\brief Tests whether the scoring scheme is set and can be invoked with the sequences passed.
    constexpr static bool expects_valid_scoring_scheme()
    {
        if constexpr (alignment_config_type::template exists<align_cfg::scoring>())
        {
            using scoring_type = std::remove_reference_t<
                                    decltype(get<align_cfg::scoring>(std::declval<alignment_config_type>()).value)
                                 >;
            return static_cast<bool>(scoring_scheme_concept<scoring_type, value_type_t<first_seq_t>,
                                                                          value_type_t<second_seq_t>>);
        }
        else
        {
            return false;
        }
    }
};

/*!\brief Configures the alignment algorithm given the sequences and the configuration object.
 * \ingroup pairwise_alignment
 */
struct alignment_configurator
{
    //!\brief Configure the algorithm.
    template <std::ranges::ViewableRange sequences_t, typename config_t>
        requires is_type_specialisation_of_v<remove_cvref_t<config_t>, configuration>
    static constexpr auto configure(sequences_t seq_range, config_t const & cfg)
    {
        // ----------------------------------------------------------------------------
        // Configure the type-erased alignment function.
        // ----------------------------------------------------------------------------

        using first_seq_t  = std::remove_reference_t<std::tuple_element_t<
                                                        0,
                                                        remove_cvref_t<decltype(*seqan3::begin(seq_range))>>>;
        using second_seq_t = std::remove_reference_t<std::tuple_element_t<
                                                        1,
                                                        remove_cvref_t<decltype(*seqan3::begin(seq_range))>>>;

        // Select the result type based on the sequences and the configuration.
        using result_t = align_result<typename align_result_selector<first_seq_t, second_seq_t, config_t>::type>;
        // Define the kernel type.
        using kernel_t = std::function<result_t(first_seq_t const &, second_seq_t const &)>;


        // ----------------------------------------------------------------------------
        // Test some basic preconditions
        // ----------------------------------------------------------------------------

        using alignment_contract_t = alignment_contract<remove_cvref_t<sequences_t>, config_t>;

        static_assert(alignment_contract_t::expects_tuple_like_value_type(),
                      "Alignment configuration error: "
                      "The value type of the sequence ranges must model the seqan3::detail::tuple_like_concept "
                      "and must contain exactly 2 elements.");

        static_assert(alignment_contract_t::expects_valid_scoring_scheme(),
                      "Alignment configuration error: "
                      "Either the scoring scheme was not configured or the given scoring scheme cannot be invoked with "
                      "the value types of the passed sequences.");

        // ----------------------------------------------------------------------------
        // Unsupported configurations
        // ----------------------------------------------------------------------------

        if constexpr (config_t::template exists<align_cfg::band>())
            throw std::domain_error{"Banded alignments are yet not supported."};

        // ----------------------------------------------------------------------------
        // Configure the algorithm
        // ----------------------------------------------------------------------------

        // Use default edit distance if gaps are not set.
        auto const & gaps = cfg.template value_or<align_cfg::gap>(gap_scheme{gap_score{-1}});
        auto const & scoring_scheme = get<align_cfg::scoring>(cfg).value;

        // Check if edit distance can be used?
        if (gaps.get_gap_open_score() == 0)
        {
            // TODO: Instead of relying on nucleotide scoring schemes we need to be able to determine the edit distance
            //       option via the scheme.
            if constexpr (is_type_specialisation_of_v<remove_cvref_t<decltype(scoring_scheme)>,
                                                      nucleotide_scoring_scheme>)
            {
                if ((scoring_scheme.score('A'_dna15, 'A'_dna15) == 0) &&
                    (scoring_scheme.score('A'_dna15, 'C'_dna15)) == -1)
                    return configure_edit_distance<kernel_t>(cfg);
            }
        }

        // ----------------------------------------------------------------------------
        // Unsupported configurations
        // ----------------------------------------------------------------------------

        if constexpr (config_t::template exists<align_cfg::aligned_ends>())
            throw std::domain_error{"Configuring the aligned ends is yet not supported."};

        if constexpr (config_t::template exists<align_cfg::result<with_begin_position_type>>())
            throw std::domain_error{"Computing the begin position is yet not supported."};

        if constexpr (config_t::template exists<align_cfg::result<with_trace_type>>())
            throw std::domain_error{"Computing the traceback is yet not supported."};

        // Configure non-edit distance alignment algorithm.
        using score_type = int;
        using cell_type = std::pair<score_type, score_type>;
        using dp_matrix_t = deferred_crtp_base<unbanded_dp_matrix_policy, std::allocator<cell_type>>;
        using affine_t = deferred_crtp_base<affine_gap_policy, cell_type>;
        using init_t = deferred_crtp_base<affine_gap_init_policy>;

        return kernel_t{alignment_algorithm<config_t, dp_matrix_t, affine_t, init_t>{cfg}};
    }

private:

    //!\brief Configure the edit distance algorithm.
    template <typename kernel_t, typename config_t>
    static constexpr auto configure_edit_distance(config_t const & cfg)
    {
        return kernel_t{edit_distance_wrapper<remove_cvref_t<config_t>>{cfg}};
    }
};

} // namespace seqan3::detail
