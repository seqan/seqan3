// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides global and local alignment configurations.
 * \author Joshua Kim <joshua.kim AT fu-berlin.de>
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 * \author Wiep van der Toorn <w.vandertoorn AT fu-berlin.de>
 * \author Lydia Buntrock <lydia.buntrock AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/detail.hpp>
#include <seqan3/core/configuration/pipeable_config_element.hpp>
#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/core/detail/strong_type.hpp>

namespace seqan3::align_cfg
{

/*!\brief Sets the local alignment method.
 * \ingroup alignment_configuration
 *
 * \details
 *
 * There are several methods for sequence alignment. We distinguish between \ref seqan3::align_cfg::method_local "local"
 * and \ref seqan3::align_cfg::method_global "global" alignments. The *semi-global* alignment is implemented as a
 * variation of the global alignment.
 *
 * \include{doc} doc/fragments/alignment_configuration_align_config_method_local.md
 *
 * ### Example
 *
 * \include test/snippet/alignment/configuration/align_cfg_method_local.cpp
 *
 * \remark For a complete overview, take a look at \ref alignment_pairwise.
 */
class method_local : private pipeable_config_element
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    method_local() = default;                                 //!< Defaulted.
    method_local(method_local const &) = default;             //!< Defaulted.
    method_local(method_local &&) = default;                  //!< Defaulted.
    method_local & operator=(method_local const &) = default; //!< Defaulted.
    method_local & operator=(method_local &&) = default;      //!< Defaulted.
    ~method_local() = default;                                //!< Defaulted.

    //!\}

    //!\privatesection
    //!\brief An internal id used to check for a valid alignment configuration.
    static constexpr seqan3::detail::align_config_id id{seqan3::detail::align_config_id::local};
};

/*!\brief A strong type representing free_end_gaps_sequence1_leading of the seqan3::align_cfg::method_global.
 * \ingroup alignment_configuration
 */
struct free_end_gaps_sequence1_leading : public seqan3::detail::strong_type<bool, free_end_gaps_sequence1_leading>
{
    //!\brief The type of the strong type base class.
    using base_t = seqan3::detail::strong_type<bool, free_end_gaps_sequence1_leading>;
    using base_t::base_t; // Import the base class constructors
};

/*!\brief A strong type representing free_end_gaps_sequence2_leading of the seqan3::align_cfg::method_global.
 * \ingroup alignment_configuration
 */
struct free_end_gaps_sequence2_leading : public seqan3::detail::strong_type<bool, free_end_gaps_sequence2_leading>
{
    //!\brief The type of the strong type base class.
    using base_t = seqan3::detail::strong_type<bool, free_end_gaps_sequence2_leading>;
    using base_t::base_t; // Import the base class constructors
};

/*!\brief A strong type representing free_end_gaps_sequence1_trailing of the seqan3::align_cfg::method_global.
 * \ingroup alignment_configuration
 */
struct free_end_gaps_sequence1_trailing : public seqan3::detail::strong_type<bool, free_end_gaps_sequence1_trailing>
{
    //!\brief The type of the strong type base class.
    using base_t = seqan3::detail::strong_type<bool, free_end_gaps_sequence1_trailing>;
    using base_t::base_t; // Import the base class constructors
};

/*!\brief A strong type representing free_end_gaps_sequence2_trailing of the seqan3::align_cfg::method_global.
 * \ingroup alignment_configuration
 */
struct free_end_gaps_sequence2_trailing : public seqan3::detail::strong_type<bool, free_end_gaps_sequence2_trailing>
{
    //!\brief The type of the strong type base class.
    using base_t = seqan3::detail::strong_type<bool, free_end_gaps_sequence2_trailing>;
    using base_t::base_t; // Import the base class constructors
};

/*!\brief Sets the global alignment method.
 * \ingroup alignment_configuration

 * \details
 *
 * There are several methods for sequence alignment. We distinguish between \ref seqan3::align_cfg::method_local "local"
 * and \ref seqan3::align_cfg::method_global "global" alignments. The *semi-global* alignment is implemented as a
 * variation of the global alignment.
 *
 * \include{doc} doc/fragments/alignment_configuration_align_config_method_global.md
 *
 * ### Example
 *
 * \include test/snippet/alignment/configuration/align_cfg_method_global.cpp
 *
 * \remark For a complete overview, take a look at \ref alignment_pairwise.
 */
class method_global : private pipeable_config_element
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    method_global() = default;                                  //!< Defaulted.
    method_global(method_global const &) = default;             //!< Defaulted.
    method_global(method_global &&) = default;                  //!< Defaulted.
    method_global & operator=(method_global const &) = default; //!< Defaulted.
    method_global & operator=(method_global &&) = default;      //!< Defaulted.
    ~method_global() = default;                                 //!< Defaulted.

    /*!\brief Construct method_global with a specific free end gap configuration.
     * \param[in] free_sequence1_leading An instance of seqan3::align_cfg::free_end_gaps_sequence1_leading that
     *                                   indicates whether leading gaps in sequence1 should be free (not penalised).
     * \param[in] free_sequence2_leading An instance of seqan3::align_cfg::free_end_gaps_sequence2_leading that
     *                                   indicates whether leading gaps in sequence2 should be free (not penalised).
     * \param[in] free_sequence1_trailing An instance of seqan3::align_cfg::free_end_gaps_sequence1_trailing that
     *                                   indicates whether trailing gaps in sequence1 should be free (not penalised).
     * \param[in] free_sequence2_trailing An instance of seqan3::align_cfg::free_end_gaps_sequence2_trailing that
     *                                   indicates whether trailing gaps in sequence2 should be free (not penalised).
     */
    constexpr method_global(seqan3::align_cfg::free_end_gaps_sequence1_leading free_sequence1_leading,
                            seqan3::align_cfg::free_end_gaps_sequence2_leading free_sequence2_leading,
                            seqan3::align_cfg::free_end_gaps_sequence1_trailing free_sequence1_trailing,
                            seqan3::align_cfg::free_end_gaps_sequence2_trailing free_sequence2_trailing) noexcept :
        free_end_gaps_sequence1_leading{free_sequence1_leading.get()},
        free_end_gaps_sequence2_leading{free_sequence2_leading.get()},
        free_end_gaps_sequence1_trailing{free_sequence1_trailing.get()},
        free_end_gaps_sequence2_trailing{free_sequence2_trailing.get()}
    {}
    //!\}

    //!\brief If set to `true`, leading gaps in sequence1 are not penalised when computing the optimal alignment.
    bool free_end_gaps_sequence1_leading{false};
    //!\brief If set to `true`, leading gaps in sequence2 are not penalised when computing the optimal alignment.
    bool free_end_gaps_sequence2_leading{false};
    //!\brief If set to `true`, trailing gaps in sequence1 are not penalised when computing the optimal alignment.
    bool free_end_gaps_sequence1_trailing{false};
    //!\brief If set to `true`, trailing gaps in sequence2 are not penalised when computing the optimal alignment.
    bool free_end_gaps_sequence2_trailing{false};

    //!\privatesection
    //!\brief An internal id used to check for a valid alignment configuration.
    static constexpr seqan3::detail::align_config_id id{seqan3::detail::align_config_id::global};
};

} // namespace seqan3::align_cfg
