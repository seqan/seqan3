// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::detail::align_config_band.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <limits>

#include <seqan3/alignment/configuration/detail.hpp>
#include <seqan3/alignment/exception.hpp>
#include <seqan3/core/configuration/pipeable_config_element.hpp>
#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/core/detail/strong_type.hpp>

namespace seqan3::align_cfg
{

/*!\brief A strong type representing the lower diagonal of the seqan3::align_cfg::band_fixed_size.
 * \ingroup alignment_configuration
 */
struct lower_diagonal : public seqan3::detail::strong_type<int32_t, lower_diagonal>
{
    //!\brief The type of the strong type base class.
    using base_t = seqan3::detail::strong_type<int32_t, lower_diagonal>;
    // Import the base class constructors
    using base_t::base_t;
};

/*!\brief A strong type representing the upper diagonal of the seqan3::align_cfg::band_fixed_size.
 * \ingroup alignment_configuration
 */
struct upper_diagonal : public seqan3::detail::strong_type<int32_t, upper_diagonal>
{
    //!\brief The type of the strong type base class.
    using base_t = seqan3::detail::strong_type<int32_t, upper_diagonal>;
    // Import the base class constructors
    using base_t::base_t;
};

/*!\brief Configuration element for setting a fixed size band.
 * \ingroup alignment_configuration
 *
 * \details
 *
 * \include{doc} doc/fragments/alignment_configuration_align_config_band.md
 *
 * ### Example
 *
 * \include test/snippet/alignment/configuration/align_cfg_band_example.cpp
 *
 * \remark For a complete overview, take a look at \ref alignment_pairwise.
 */
class band_fixed_size : private pipeable_config_element
{
public:
    //!\brief The selected lower diagonal. Defaults to `std::%numeric_limits<int32_t>::%lowest()`.
    int32_t lower_diagonal{std::numeric_limits<int32_t>::lowest()};
    //!\brief The selected upper diagonal. Defaults to `std::%numeric_limits<int32_t>::%max()`.
    int32_t upper_diagonal{std::numeric_limits<int32_t>::max()};

    /*!\name Constructor, destructor and assignment
     * \{
     */
    constexpr band_fixed_size() = default;                                    //!< Defaulted.
    constexpr band_fixed_size(band_fixed_size const &) = default;             //!< Defaulted.
    constexpr band_fixed_size(band_fixed_size &&) = default;                  //!< Defaulted.
    constexpr band_fixed_size & operator=(band_fixed_size const &) = default; //!< Defaulted.
    constexpr band_fixed_size & operator=(band_fixed_size &&) = default;      //!< Defaulted.
    ~band_fixed_size() = default;                                             //!< Defaulted.

    /*!\brief Initialises the fixed size band by setting the lower and the upper matrix diagonal.
     *
     * \param lower_diagonal \copybrief seqan3::align_cfg::band_fixed_size::lower_diagonal
     * \param upper_diagonal \copybrief seqan3::align_cfg::band_fixed_size::upper_diagonal
     *
     * \details
     *
     * The lower diagonal represents the lower bound of the banded matrix, i.e. the alignment cannot pass below this
     * diagonal. Similar, the upper diagonal represents the upper bound of the alignment. During the alignment
     * configuration and execution the band parameters will be checked and an exception will be thrown in case of
     * an invalid configuration.
     */
    constexpr band_fixed_size(seqan3::align_cfg::lower_diagonal const lower_diagonal,
                              seqan3::align_cfg::upper_diagonal const upper_diagonal) :
        lower_diagonal{lower_diagonal.get()},
        upper_diagonal{upper_diagonal.get()}
    {}
    //!\}

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr seqan3::detail::align_config_id id{seqan3::detail::align_config_id::band};
};

} // namespace seqan3::align_cfg
