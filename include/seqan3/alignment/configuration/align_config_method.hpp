// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

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
 * and \ref seqan3::align_cfg::method_global "global" alignments and the semi-global alignment in between, which is a
 * variation of the global alignment.
 *
 * **Global Alignment**:\verbatim
--T--CC-C-AGT--TATGT-CAGGGGACACG-A-GCATGCAGA-GAC
  |  || |  ||  | | | |||    || | | |  | ||||   |
AATTGCCGCC-GTCGT-T-TTCAG----CA-GTTATG-T-CAGAT--C\endverbatim
 * Finding the optimal global alignment of two sequences is solved by the **Needleman-Wunsch algorithm** a “dynamic
 * program”. Hence, the algorithm requires O(nm) time and memory.
 *
 * **Local Alignment** (better suited to find conserved segment):
 * \verbatim
                  tccCAGTTATGTCAGgggacacgagcatgcagagac
                     ||||||||||||
aattgccgccgtcgttttcagCAGTTATGTCAGatc
\endverbatim
 * A \ref seqan3::align_cfg::method_local "local" alignment is effectively a global alignment of two partial sequences,
 * that is, when two genes from different species are similar in short conserved regions and dissimilar in the remaining
 * regions. A global alignment would not find the local matching because it would try to align the entire sequence.
 * This is solved by the **Smith-Waterman algorithm**.
 *
 * **Complexity** of both algorithms:
 *
 * We need to store (n+1)(m+1) numbers, with a constant number of calculations required for each number (three sums and
 * one maximum). Therefore, the algorithm requires **O(nm)** time and memory.
 *
 * To reduce the time complexity you can use a \ref seqan3::align_cfg::band_fixed_size "banded" alignment.
 *
 * --
 *
 * By distinguishing between these alignments, it is possible to define a subset of
 * \ref seqan3::align_cfg "configurations" that can work with a particular method. Since it is not possible to guess
 * which method a user wants, there is no default and this \ref seqan3::align_cfg "configuration" must always be
 * specified for the alignment algorithm.
 *
 * These include configurations for a semi-global alignment, such as the overlap alignment when you do not wish to
 * penalize \ref seqan3::align_cfg::free_end_gaps_sequence1_leading "first sequence leading",
 * \ref seqan3::align_cfg::free_end_gaps_sequence2_leading "second sequence leading" or
 * \ref seqan3::align_cfg::free_end_gaps_sequence1_trailing "first sequence trailing",
 * \ref seqan3::align_cfg::free_end_gaps_sequence2_trailing "second sequence trailing" gaps.
 *
 * If you want to define affine gap costs instead of linear ones, you can define a
 * \ref seqan3::align_cfg::gap_cost_affine "affine gap cost scheme". If the gap scheme is not configured, it will
 * default to a linear gap scheme initialised with edit distance.
 *
 * \see [lecture script - pairwise alignment](https://www.mi.fu-berlin.de/en/inf/groups/abi/teaching/lectures/lectures_past/WS0910/V___Algorithmen_und_Datenstrukturen/scripts/alignment.pdf)
 *
 * ### Example
 *
 * \include test/snippet/alignment/configuration/align_cfg_method_local.cpp
 */
class method_local : private pipeable_config_element
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    method_local() = default; //!< Defaulted.
    method_local(method_local const &) = default; //!< Defaulted.
    method_local(method_local &&) = default; //!< Defaulted.
    method_local & operator=(method_local const &) = default; //!< Defaulted.
    method_local & operator=(method_local &&) = default; //!< Defaulted.
    ~method_local() = default; //!< Defaulted.

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
 * \copydetails seqan3::align_cfg::method_local
 */
class method_global : private pipeable_config_element
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    method_global() = default; //!< Defaulted.
    method_global(method_global const &) = default; //!< Defaulted.
    method_global(method_global &&) = default; //!< Defaulted.
    method_global & operator=(method_global const &) = default; //!< Defaulted.
    method_global & operator=(method_global &&) = default; //!< Defaulted.
    ~method_global() = default; //!< Defaulted.

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
