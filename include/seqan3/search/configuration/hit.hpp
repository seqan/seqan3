// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides the configuration to define the hit strategies "hit_strata", "hit_all",
 *        "hit_all_best", "hit_single_best".
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <variant>

#include <seqan3/core/configuration/pipeable_config_element.hpp>
#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/search/configuration/detail.hpp>
#include <seqan3/utility/type_pack/traits.hpp>

namespace seqan3::search_cfg
{

/*!\brief Configuration element to receive all hits within the error bounds.
 * \ingroup search_configuration
 * \see search_configuration
 * \sa \ref search_configuration_subsection_hit_strategy "Section on Hit Strategy"
 */
class hit_all : private pipeable_config_element
{
public:
    /*!\name Constructors, assignment and destructor
     * \{
     */
    constexpr hit_all() = default;                            //!< Defaulted.
    constexpr hit_all(hit_all const &) = default;             //!< Defaulted.
    constexpr hit_all(hit_all &&) = default;                  //!< Defaulted.
    constexpr hit_all & operator=(hit_all const &) = default; //!< Defaulted.
    constexpr hit_all & operator=(hit_all &&) = default;      //!< Defaulted.
    ~hit_all() = default;                                     //!< Defaulted.

    //!\}

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::search_config_id id{detail::search_config_id::hit};
};

/*!\brief Configuration element to receive all hits with the lowest number of errors within the error bounds.
 * \ingroup search_configuration
 * \see search_configuration
 * \sa \ref search_configuration_subsection_hit_strategy "Section on Hit Strategy"
 */
class hit_all_best : private pipeable_config_element
{
public:
    /*!\name Constructors, assignment and destructor
     * \{
     */
    constexpr hit_all_best() = default;                                 //!< Defaulted.
    constexpr hit_all_best(hit_all_best const &) = default;             //!< Defaulted.
    constexpr hit_all_best(hit_all_best &&) = default;                  //!< Defaulted.
    constexpr hit_all_best & operator=(hit_all_best const &) = default; //!< Defaulted.
    constexpr hit_all_best & operator=(hit_all_best &&) = default;      //!< Defaulted.
    ~hit_all_best() = default;                                          //!< Defaulted.

    //!\}

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr seqan3::detail::search_config_id id{seqan3::detail::search_config_id::hit};
};

/*!\brief Configuration element to receive a single best hit with the lowest number of errors within the error bounds.
 * \ingroup search_configuration
 * \see search_configuration
 * \sa \ref search_configuration_subsection_hit_strategy "Section on Hit Strategy"
 */
class hit_single_best : private pipeable_config_element
{
public:
    /*!\name Constructors, assignment and destructor
     * \{
     */
    constexpr hit_single_best() = default;                                    //!< Defaulted.
    constexpr hit_single_best(hit_single_best const &) = default;             //!< Defaulted.
    constexpr hit_single_best(hit_single_best &&) = default;                  //!< Defaulted.
    constexpr hit_single_best & operator=(hit_single_best const &) = default; //!< Defaulted.
    constexpr hit_single_best & operator=(hit_single_best &&) = default;      //!< Defaulted.
    ~hit_single_best() = default;                                             //!< Defaulted.

    //!\}

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr seqan3::detail::search_config_id id{seqan3::detail::search_config_id::hit};
};

/*!\brief Configuration element to receive all hits with the best number of errors plus the given stratum.
 *        All hits are found with the fewest number of errors plus 'stratum'.
 * \ingroup search_configuration
 * \see search_configuration
 * \sa \ref search_configuration_subsection_hit_strategy "Section on Hit Strategy"
 */
class hit_strata : private pipeable_config_element
{
public:
    //!\brief The stratum value [default: 0].
    uint8_t stratum{};

    /*!\name Constructors, assignment and destructor
     * \{
     */
    constexpr hit_strata() = default;                               //!< Defaulted.
    constexpr hit_strata(hit_strata const &) = default;             //!< Defaulted.
    constexpr hit_strata(hit_strata &&) = default;                  //!< Defaulted.
    constexpr hit_strata & operator=(hit_strata const &) = default; //!< Defaulted.
    constexpr hit_strata & operator=(hit_strata &&) = default;      //!< Defaulted.
    ~hit_strata() = default;                                        //!< Defaulted.

    /*!\brief Initialises the strata config.
     * \param[in] stratum The stratum to include in the search.
     */
    hit_strata(uint32_t stratum) : stratum{static_cast<uint8_t>(stratum)}
    {}
    //!\}

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::search_config_id id{detail::search_config_id::hit};
};

/*!\brief A dynamic configuration element to configure the hit strategy at runtime.
 * \ingroup search_configuration
 * \see search_configuration
 * \sa \ref search_configuration_subsection_hit_strategy "Section on Hit Strategy"
 */
class hit : private pipeable_config_element
{
public:
    /*!\brief The type of the std::variant holding the hit configuration element alternatives
     *
     * The additional detail::empty_type marks the seqan3::search_cfg::hit as default constructed, with no selected hit
     * configuration which can be checked within the search algorithm.
     */
    using hit_variant_type = std::variant<detail::empty_type, hit_all, hit_all_best, hit_single_best, hit_strata>;

    /*!\name Constructors, assignment and destructor
     * \{
     */
    hit() = default;                        //!< Defaulted.
    hit(hit const &) = default;             //!< Defaulted.
    hit(hit &&) = default;                  //!< Defaulted.
    hit & operator=(hit const &) = default; //!< Defaulted.
    hit & operator=(hit &&) = default;      //!< Defaulted.
    ~hit() = default;                       //!< Defaulted.

    /*!\brief Sets the given configuration element to the dynamic hit configuration element.
     *
     * \tparam hit_config_t The type of the hit configuration element to set.
     * \param[in] hit_config The hit configuration element to set.
     *
     * \details
     *
     * Only one of the static hit configuration elements are valid: seqan3::search_cfg::hit_all,
     * seqan3::search_cfg::hit_all_best, seqan3::search_cfg::hit_single_best and seqan3::search_cfg::hit_strata.
     */
    template <typename hit_config_t>
        requires pack_traits::contains<hit_config_t, hit_all, hit_all_best, hit_single_best, hit_strata>
    explicit hit(hit_config_t hit_config) noexcept : hit_variant{std::move(hit_config)}
    {}

    //!\copydoc seqan3::search_cfg::hit::hit(hit_config_t hit_config)
    template <typename hit_config_t>
        requires pack_traits::contains<hit_config_t, hit_all, hit_all_best, hit_single_best, hit_strata>
    hit & operator=(hit_config_t hit_config) noexcept
    {
        hit_variant = std::move(hit_config);
        return *this;
    }
    //!\}

    //!\brief A std::variant over the valid hit configuration elements.
    hit_variant_type hit_variant{};

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::search_config_id id{detail::search_config_id::hit};
};

} // namespace seqan3::search_cfg
