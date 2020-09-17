// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides some detail functionality to the algorithm seqan3::align_multiple.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \author Joerg Winkler <j.winkler AT fu-berlin.de>
 * \author Simon Sasse <simon.sasse AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/algorithm>
#include <cstdlib>
#include <string>
#include <vector>

#include <seqan3/core/platform.hpp>

#if SEQAN3_HAS_SEQAN2
#include <seqan/graph_msa.h>
#include <seqan/reduced_aminoacid.h>
#endif

#include <seqan3/alignment/configuration/align_config_band.hpp>
#include <seqan3/alignment/configuration/align_config_gap.hpp>
#include <seqan3/alignment/configuration/align_config_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/aminoacid/all.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/core/detail/test_accessor.hpp>
#include <seqan3/core/type_list/traits.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/range/views/char_to.hpp>
#include <seqan3/range/views/chunk.hpp>

namespace seqan3::detail
{

/*!\brief Holds all functionality to make the seqan3::align_multiple interface compatible with the SeqAn2 msa algorithm.
 * \ingroup multiple_alignment
 * \sa seqan3::align_multiple
 */
template <typename seqan3_alphabet_type>
class align_multiple_seqan2_adaptation
{
private:
    //!\brief A list of SeqAn3 types that are currently supported for seqan3::align_multiple.
    using seqan3_types = type_list<dna4, dna5, dna15, rna4, rna5, aa27, aa10murphy, aa10li>;
    //!\brief The corresponding list of SeqAn2 types (the order must be the same!).
    using seqan2_types = type_list<seqan::Dna,
                                   seqan::Dna5,
                                   seqan::Iupac,
                                   seqan::Rna,
                                   seqan::Rna5,
                                   seqan::AminoAcid,
                                   seqan::ReducedAminoAcid<seqan::Murphy10>,
                                   seqan::ReducedAminoAcid<seqan::Li10>>;
    //!\brief The index of the seqan3_alphabet_type in the list of seqan3_types.
    static constexpr auto index = list_traits::find<seqan3_alphabet_type, seqan3_types>;

    static_assert(index != -1, "The seqan3_alphabet_type is currently not supported.");

public:
    //!\brief The (SeqAn2) alphabet type used in the multiple sequence alignment algorithm.
    using alphabet_type = list_traits::at<index, seqan2_types>;
    //!\brief The (SeqAn2) sequence tupe used in the multiple sequence alignment algorithm.
    using sequence_type = seqan::String<alphabet_type>;
    //!\brief The output graph type of the multiple sequence alignment algorithm.
    using graph_type = seqan::Graph<seqan::Alignment<seqan::StringSet<sequence_type, seqan::Dependent<>>,
                                                     void,
                                                     seqan::WithoutEdgeId>>;

    /*!\brief Create the SeqAn2 configuration object based on the given SeqAn3 configuration.
     * \tparam seqan3_configuration_t The type of the SeqAn3 configuration object.
     * \param[in] config The SeqAn3 configuration object.
     * \return A SeqAn2 *MsaOptions* object equivalent to the SeqAn3 configuration.
     */
    template <typename seqan3_configuration_t>
    auto create_msa_configuration(seqan3_configuration_t const & config)
    {
        validate_configuration(config);

        auto msa_options = initialise_scoring_scheme(config);

        seqan::appendValue(msa_options.method, 0); // global alignment
        seqan::appendValue(msa_options.method, 1); // local alignment
        msa_options.build = 0; // neighbour joining to build the guide tree

        if constexpr (config.template exists<seqan3::align_cfg::band_fixed_size>())
        {
            msa_options.pairwiseAlignmentMethod = 2; // banded
            auto const & band = get<seqan3::align_cfg::band_fixed_size>(config);
            msa_options.bandWidth = band.upper_diagonal.get() - band.lower_diagonal.get();
        }
        else
        {
            msa_options.pairwiseAlignmentMethod = 1; // unbanded
        }

        // SeqAn2 tcoffee app default: gap -1, gap open -13
        auto const & gaps = config.get_or(align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-13}}}).value;
        // convert to SeqAn2 gap score convention.
        // See: https://docs.seqan.de/seqan/3-master-user/classseqan3_1_1gap__scheme.html
        msa_options.sc.data_gap_open = gaps.get_gap_open_score() + gaps.get_gap_score();
        msa_options.sc.data_gap_extend = gaps.get_gap_score();

        return msa_options;
    }

    /*!\brief Convert a given range of SeqAn3 sequences to a SeqAn2 StringSet of sequences.
     * \tparam range_t The type of SeqAn3 range, must model std::ranges::forward_range.
     * \param input The range of SeqAn3 sequences.
     * \return A std::pair of (1) a SeqAn2 *StringSet* container holding the sequences and (2) another
     *  SeqAn2 *StringSet* container holding the sequence ids.
     */
    template <std::ranges::forward_range range_t>
    auto convert_sequences(std::vector<range_t> const & input)
    {
        seqan::StringSet<sequence_type, seqan::Owner<>> sequences;
        seqan::StringSet<seqan::String<char>> ids;

        seqan::String<char> const dummy_name{"dummy_name"};
        sequence_type tmp;

        for (auto const & seq : input)
        {
            for (auto chr : seq)
                seqan::appendValue(tmp, alphabet_type{seqan3::to_char(chr)});

            seqan::appendValue(sequences, tmp);
            seqan::appendValue(ids, dummy_name);
            seqan::clear(tmp);
        }

        return std::make_pair(std::move(sequences), std::move(ids));
    }

    /*!\brief Create a vector of gapped SeqAn3 sequences from a SeqAn2 alignment graph.
     * \tparam graph_type The type of the alignment graph resulting from a seqan::globalMsaAlignment.
     * \param alignment_graph The alignment graph from a seqan::globalMsaAlignment.
     * \return A std::vector of gapped SeqAn3 sequences representing the MSA result.
     */
    template <typename graph_type>
    auto create_output(graph_type & alignment_graph)
    {
        using gapped_alphabet_type = seqan3::gapped<seqan3_alphabet_type>;

        // convert the alignment graph into a string of all aligned sequences concatenated to each other
        std::string concatenated_alignment_string;
        seqan::convertAlignment(alignment_graph, concatenated_alignment_string);

        // initialise the output vector of SeqAn3 sequences
        auto & underlying_string_set = seqan::stringSet(alignment_graph);
        std::vector<std::vector<gapped_alphabet_type>> output{seqan::length(underlying_string_set)};

        size_t const alignment_length = concatenated_alignment_string.size() / output.size();

        auto iter = output.begin();
        for (auto && gapped_seq : concatenated_alignment_string | seqan3::views::char_to<gapped_alphabet_type>
                                                                | seqan3::views::chunk(alignment_length))
        {
            iter->reserve(alignment_length);
            std::ranges::copy(gapped_seq, std::cpp20::back_inserter(*iter++));
        }

        return output;
    }

private:
    /*!\brief Predicate that represents whether `candidate_t` is a type that is not allowed as MSA configuration.
     * \tparam candidate_t The type to test.
     *
     * \details
     * The predicate's value is set to: NOT (band OR gap OR scoring).
     */
    template <typename candidate_t>
    struct is_invalid_msa_config :
        std::negation<std::disjunction<std::is_same<candidate_t, seqan3::align_cfg::band_fixed_size>,
                                       is_type_specialisation_of<candidate_t, seqan3::align_cfg::gap>,
                                       is_type_specialisation_of<candidate_t, seqan3::align_cfg::scoring_scheme>>>{};

    /*!\brief Validate the given MSA configuration.
     * \tparam configs_t The specified configuration elements.
     */
    template <typename ... configs_t>
    void validate_configuration(seqan3::configuration<configs_t...> const &)
    {
        static_assert(seqan3::pack_traits::find_if<is_invalid_msa_config, configs_t...> == -1,
                      "The given MSA configuration is not valid.");
    }

    /*!\brief Create the SeqAn2 scoring scheme based on the given SeqAn3 configuration.
     * \tparam seqan3_configuration_t The type of the configuration object.
     * \param config The configuration that contains scoring parameters.
     * \return A SeqAn2 *MsaOptions* object with initialised scores.
     *
     * \details
     * The scoring scheme matrices are copied from SeqAn3 into the corresponding SeqAn2 type.
     * Note that the alphabet type can be different to that of the input sequences
     * (seqan3::detail::align_multiple_seqan2_adaptation::alphabet_type), e.g., the SeqAn3 nucleotide scoring scheme is
     * defined for dna16 but works for most other nucleotide alphabets too. The same holds for the SeqAn2 equivalent.
     */
    template <typename seqan3_configuration_t>
    //!\cond
        requires (seqan3_configuration_t::template exists<seqan3::align_cfg::scoring_scheme>())
    //!\endcond
    auto initialise_scoring_scheme(seqan3_configuration_t const & config)
    {
        auto scoring_scheme = get<seqan3::align_cfg::scoring_scheme>(config).value;

        // A SeqAn3 scoring scheme might be of alphabet type dna15 although the aligned alphabet is dna4
        // SeqAn3 alphabet type used:
        using scoring_scheme_alphabet_type = typename decltype(scoring_scheme)::alphabet_type;
        constexpr auto scoring_scheme_alphabet_index = list_traits::find<scoring_scheme_alphabet_type, seqan3_types>;
        // SeqAn2 alphabet type to use accordingly:
        using score_matrix_alphabet_type = list_traits::at<scoring_scheme_alphabet_index, seqan2_types>;
        using score_matrix_type = seqan::Score<int32_t, seqan::ScoreMatrix<score_matrix_alphabet_type>>;

        seqan::MsaOptions<alphabet_type, score_matrix_type> msa_options{};
        score_matrix_type score_matrix;

        for (size_t i = 0; i < seqan3::alphabet_size<scoring_scheme_alphabet_type>; ++i)
        {
            for (size_t j = 0; j < seqan3::alphabet_size<scoring_scheme_alphabet_type>; ++j)
            {
                auto seqan3_i = seqan3::assign_rank_to(i, scoring_scheme_alphabet_type{});
                auto seqan3_j = seqan3::assign_rank_to(j, scoring_scheme_alphabet_type{});
                score_matrix_alphabet_type seqan2_i{seqan3::to_char(seqan3_i)};
                score_matrix_alphabet_type seqan2_j{seqan3::to_char(seqan3_j)};

                setScore(score_matrix, seqan2_i, seqan2_j, scoring_scheme.score(seqan3_i, seqan3_j));
            }
        }

        msa_options.sc = score_matrix;
        return msa_options;
    }

    /*!\brief Create the SeqAn2 scoring scheme based on default values.
     * \tparam seqan3_configuration_t The type of a configuration object which does not contain scoring parameters.
     * \return A SeqAn2 *MsaOptions* object with initialised scores.
     *
     * Settings:
     *   * scoring for amino acid sequences: Blosum62 matrix
     *   * scoring for nucleotide sequences: match = +5 and mismatch = -4
     */
    template <typename seqan3_configuration_t>
    //!\cond
        requires (!seqan3_configuration_t::template exists<seqan3::align_cfg::scoring_scheme>())
    //!\endcond
    auto initialise_scoring_scheme(seqan3_configuration_t const &)
    {
        using score_type = std::conditional_t<std::same_as<alphabet_type, seqan::AminoAcid> ||
                                              std::same_as<alphabet_type, seqan::ReducedAminoAcid<seqan::Murphy10>> ||
                                              std::same_as<alphabet_type, seqan::ReducedAminoAcid<seqan::Li10>>,
                                              seqan::Blosum62,
                                              seqan::Score<int>>;

        seqan::MsaOptions<alphabet_type, score_type> msa_options{};

        if constexpr (std::is_same_v<score_type, seqan::Score<int>>)
        {
            msa_options.sc.data_match = 5;       // tcoffee app default
            msa_options.sc.data_mismatch = -4;   // tcoffee app default
        }

        return msa_options;
    }

    //!\brief Befriend seqan3::detail::test_accessor to grant access to private member functions.
    friend struct ::seqan3::detail::test_accessor;
};

} // namespace seqan3::detail
