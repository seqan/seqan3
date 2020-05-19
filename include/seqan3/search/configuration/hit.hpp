// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the configuration to define the hit strategies "hit_strata", "hit_all",
 *        "hit_all_best", "hit_single_best".
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <variant>

#include <seqan3/core/algorithm/pipeable_config_element.hpp>
#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/search/configuration/detail.hpp>

namespace seqan3::detail
{

/*!\brief Configuration element to receive all hits within the error bounds.
 * \ingroup search_configuration
 */
struct hit_all_tag : public pipeable_config_element<hit_all_tag, empty_type>
{
    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::search_config_id id{detail::search_config_id::hit};
};

/*!\brief Configuration element to receive all hits with the lowest number of errors within the error bounds.
 * \ingroup search_configuration
 */
struct hit_all_best_tag : public pipeable_config_element<hit_all_best_tag, empty_type>
{
    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::search_config_id id{detail::search_config_id::hit};
};

/*!\brief Configuration element to receive a single best hit with the lowest number of errors within the error bounds.
 * \ingroup search_configuration
 */
struct hit_single_best_tag : public pipeable_config_element<hit_single_best_tag, empty_type>
{
    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::search_config_id id{detail::search_config_id::hit};
};

} // namespace seqan3::detail

namespace seqan3::search_cfg
{

/*!\brief \copybrief seqan3::detail::hit_all_tag
 * \ingroup search_configuration
 * \sa \ref search_configuration_subsection_hit_strategy "Section on Hit Strategy"
 * \hideinitializer
 */
inline constexpr detail::hit_all_tag hit_all{};

/*!\brief \copybrief seqan3::detail::hit_all_best_tag
 * \ingroup search_configuration
 * \sa \ref search_configuration_subsection_hit_strategy "Section on Hit Strategy"
 * \hideinitializer
 */
inline constexpr detail::hit_all_best_tag hit_all_best{};

/*!\brief \copybrief seqan3::detail::hit_single_best_tag
 * \ingroup search_configuration
 * \sa \ref search_configuration_subsection_hit_strategy "Section on Hit Strategy"
 * \hideinitializer
 */
inline constexpr detail::hit_single_best_tag hit_single_best{};

/*!\brief Configuration element to receive all hits with the best number of errors plus the given stratum.
 *        All hits are found with the fewest number of errors plus 'stratum'.
 * \ingroup search_configuration
 * \sa \ref search_configuration_subsection_hit_strategy "Section on Hit Strategy"
 */
struct hit_strata : public pipeable_config_element<hit_strata, uint8_t>
{
    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::search_config_id id{detail::search_config_id::hit};
};

/*!\brief A dynamic configuration element to select the search strategy.
 * \ingroup search_configuration
 *
 * todo
 *
 * \sa \ref search_configuration_subsection_hit_strategy "Section on Hit Strategy"
 */
struct hit : public pipeable_config_element<hit, std::variant<seqan3::detail::hit_all_tag,
                                                              seqan3::detail::hit_all_best_tag,
                                                              seqan3::detail::hit_single_best_tag,
                                                              hit_strata>>
{
private:
    //!\brief A variant including all possible search hit types.
    using variant_t = std::variant<seqan3::detail::hit_all_tag,
                                   seqan3::detail::hit_all_best_tag,
                                   seqan3::detail::hit_single_best_tag,
                                   hit_strata>;
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    hit() = default; //!< Defaulted.
    hit(hit const &) = default; //!< Defaulted.
    hit & operator=(hit const &) = default; //!< Defaulted.
    hit(hit &&) = default; //!< Defaulted.
    hit & operator=(hit &&) = default; //!< Defaulted.
    ~hit() = default; //!< Defaulted.

    /*!\brief Construct from a search strategy, e.g. seqan3::search_cfg::hit_strata.
     * \tparam option_t The type of search strategy.
     * \param[in] option The search strategy to select.
     */
    template <typename option_t>
    //!\cond
        requires seqan3::list_traits::contains<option_t, detail::transfer_template_args_onto_t<variant_t, type_list>>
    //!\endcond
    hit(option_t option) : selection{std::move(option)}
    {}

    /*!\brief Assign from a search strategy, e.g. seqan3::search_cfg::hit_strata.
     * \tparam option_t The type of search strategy.
     * \param[in] option The search strategy to select.
     */
    template <typename option_t>
    //!\cond
        requires seqan3::list_traits::contains<option_t, detail::transfer_template_args_onto_t<variant_t, type_list>>
    //!\endcond
    hit & operator=(option_t option)
    {
        selection = std::move(option);
        return *this;
    }
    //!\}

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::search_config_id id{detail::search_config_id::hit};

    //!\brief The selected search strategy.
    variant_t selection{};
};

} // namespace seqan3::search_cfg
