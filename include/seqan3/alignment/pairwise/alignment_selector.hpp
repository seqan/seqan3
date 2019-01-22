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
#include <seqan3/alignment/pairwise/alignment_algorithm.hpp>
#include <seqan3/alignment/pairwise/align_result.hpp>
#include <seqan3/alignment/pairwise/edit_distance_unbanded.hpp>
#include <seqan3/alignment/pairwise/policy/affine_gap_init_policy.hpp>
#include <seqan3/alignment/pairwise/policy/affine_gap_policy.hpp>
#include <seqan3/alignment/pairwise/policy/unbanded_dp_matrix_policy.hpp>
#include <seqan3/alignment/pairwise/align_result_selector.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/core/metafunction/deferred_crtp_base.hpp>
#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/core/metafunction/template_inspection.hpp>
#include <seqan3/core/type_list.hpp>
#include <seqan3/range/view/persist.hpp>

namespace seqan3::detail
{

/*!\brief Configures the alignment kernel given the sequences and the configuration object.
 * \ingroup pairwise_alignment
 */
struct alignment_configurator
{
    //!\brief Configure the edit distance algorithm.
    template <typename kernel_t, typename config_t>
    static constexpr auto configure_edit_distance(config_t const & cfg)
    {
        return kernel_t{edit_distance_wrapper<remove_cvref_t<config_t>>{cfg}};
    }

    //!\brief Configure the algorithm.
    template <std::ranges::View sequences_t, typename config_t>
        requires is_type_specialisation_of_v<remove_cvref_t<config_t>, configuration>
    static constexpr auto configure(sequences_t SEQAN3_DOXYGEN_ONLY(seq_range), config_t const & cfg)
    {
        using first_seq_t = std::remove_reference_t<
                                std::tuple_element_t<
                                    0,
                                    value_type_t<std::ranges::iterator_t<remove_cvref_t<sequences_t>>>
                                >
                            >;
        using second_seq_t = std::remove_reference_t<
                                std::tuple_element_t<
                                    1,
                                    value_type_t<std::ranges::iterator_t<remove_cvref_t<sequences_t>>>
                                >
                             >;

        using result_t = align_result<typename align_result_selector<first_seq_t,
                                                                     second_seq_t,
                                                                     remove_cvref_t<config_t>>::type
                                     >;
        using kernel_t = std::function<result_t(first_seq_t const &, second_seq_t const &)>;

        auto const & gaps = cfg.template value_or<align_cfg::gap>(gap_scheme{gap_score{-1}});
        auto const & scoring_scheme =
            cfg.template value_or<align_cfg::scoring>(nucleotide_scoring_scheme{match_score{0}, mismatch_score{-1}});
        // Linear gaps
        if (gaps.get_gap_open_score() == 0)
        {
            if constexpr (is_type_specialisation_of_v<remove_cvref_t<decltype(scoring_scheme)>,
                                                      nucleotide_scoring_scheme>)
            {
                // TODO: Check if the matrix is Levenshtein distance.
                if ((scoring_scheme.score('A'_dna15, 'A'_dna15) == 0) &&
                    (scoring_scheme.score('A'_dna15, 'C'_dna15)) == -1)
                    return configure_edit_distance<kernel_t>(cfg);
                else
                    throw std::domain_error{"Linear gaps are not yet implemented."};
            }
            else // Not nucleotide_scoring_scheme
            {
                throw std::domain_error{"Linear gaps are not yet implemented."};
            }
        }
        else // Affine gaps
        {
            using score_type = int;
            using cell_type = std::pair<score_type, score_type>;
            using dp_matrix_t = deferred_crtp_base<unbanded_dp_matrix_policy, std::allocator<cell_type>>;
            using affine_t = deferred_crtp_base<affine_gap_policy, cell_type>;
            using init_t = deferred_crtp_base<affine_gap_init_policy>;
            //TODO: copies but we want to have unique_funtion to only move.
            return kernel_t{alignment_algorithm<config_t, dp_matrix_t, affine_t, init_t>{cfg}};
        }
    }
};

} // namespace seqan3::detail
