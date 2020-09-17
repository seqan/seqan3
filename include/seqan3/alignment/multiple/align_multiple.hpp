// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the algorithm seqan3::align_multiple for multiple sequence alignment.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \author Joerg Winkler <j.winkler AT fu-berlin.de>
 * \author Simon Sasse <simon.sasse AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/ranges>
#include <vector>

#include <seqan3/core/platform.hpp>

#if SEQAN3_HAS_SEQAN2
#include <seqan/graph_msa.h>
#endif

#include <seqan3/alignment/configuration/align_config_gap.hpp>
#include <seqan3/alignment/multiple/detail/align_multiple_seqan2_adaptation.hpp>

namespace seqan3::align_cfg
{

/*!\brief The standard configuration for multiple sequence alignment.
 * \ingroup multiple_alignment
 * \details
 *
 * The standard configuration values are adapted from SeqAn2 in order to provide the same behaviour.
 *
 *   * general gap score: -1
 *   * additional gap open score: -13
 *   * band constraints: none
 *   * scoring for amino acid sequences: Blosum62 matrix (set internally)
 *   * scoring for nucleotide sequences: match = +5 and mismatch = -4 (set internally)
 *
 * \if DEV
 * See seqan3::detail::align_multiple_seqan2_adaptation::create_msa_configuration for more details on the configuration.
 * \endif
 */
constexpr configuration msa_default_configuration = gap{gap_scheme{gap_score{-1}, gap_open_score{-13}}};

} // namespace seqan3::align_cfg

namespace seqan3
{

/*!\brief The algorithm for multiple sequence alignment.
 * \ingroup multiple_alignment
 * \tparam range_t Type of the input sequences, must model std::ranges::forward_range.
 * \tparam config_t Type of the configuration; defaults to the type of seqan3::align_cfg::msa_default_configuration.
 * \param[in] input A vector of sequences that you want to align.
 * \param[in] config A configuration object that stores the settings for the algorithm;
 *               defaults to seqan3::align_cfg::msa_default_configuration.
 * \return The multiple sequence alignment as a vector of gapped sequences.
 *
 * \details
 *
 * Computes a multiple sequence alignment from the given input sequences, using a consistency-based progressive
 * alignment algorithm on a graph of sequence segments. You can use the configuration object to
 * specify various parameters, like gap scores, alignment scores and band constraints. The return type is
 * `std::vector<std::vector<gapped_alphabet_type>>`, with the inner type derived from the input sequence type.
 */
template <std::ranges::forward_range range_t, typename config_t = decltype(align_cfg::msa_default_configuration)>
auto align_multiple(std::vector<range_t> const & input, config_t config = align_cfg::msa_default_configuration)
{
#if !SEQAN3_HAS_SEQAN2 // multiple sequence alignment is only enabled with seqan2
    static_assert(false, "You need to have at least seqan 2.4");
#else // SEQAN3_HAS_SEQAN2

    using seqan3_alphabet_type = std::ranges::range_value_t<range_t>;
    using seqan2_adaptation_type = detail::align_multiple_seqan2_adaptation<seqan3_alphabet_type>;
    using seqan2_graph_type = typename seqan2_adaptation_type::graph_type;

    seqan2_adaptation_type seqan2_adaptation{};

    auto seqan2_msa_options = seqan2_adaptation.create_msa_configuration(config);
    auto && [sequences, ids] = seqan2_adaptation.convert_sequences(input);

    seqan2_graph_type alignment_graph;

    seqan::globalMsaAlignment(alignment_graph, sequences, ids, seqan2_msa_options);

    return seqan2_adaptation.create_output(alignment_graph);

#endif // !SEQAN3_HAS_SEQAN2
}

} // namespace seqan3
