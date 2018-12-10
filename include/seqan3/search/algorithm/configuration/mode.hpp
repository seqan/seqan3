// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

/*!\file
 * \brief Provides the mode configuration to define the search modes "all", "all_best", "best" and "strata".
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/algorithm/pipeable_config_element.hpp>
#include <seqan3/core/detail/strong_type.hpp>
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/search/algorithm/configuration/detail.hpp>

/*!\addtogroup search
 * \{
 */

namespace seqan3::detail
{

//!\brief Type for the "all" value for the configuration element "mode".
struct search_mode_all {};
//!\brief Type for the "all_best" value for the configuration element "mode".
struct search_mode_all_best {};
//!\brief Type for the "best" value for the configuration element "mode".
struct search_mode_best {};

} // namespace seqan3::detail

namespace seqan3::search_cfg
{

//!\brief Configuration element to receive all hits within the error bounds.
inline detail::search_mode_all constexpr all;
//!\brief Configuration element to receive all hits within the lowest number of errors.
inline detail::search_mode_all_best constexpr all_best;
//!\brief Configuration element to receive one best hit (with the lowest number of errors).
inline detail::search_mode_best constexpr best;

/*!\brief Configuration element to receive all hits with the best number of errors plus the strata value.
 *        A strong type of underlying type `uint8_t` that represents the number or errors for strata.
 *        All hits are found with the fewest numbererrors plus 'value'.
 * \tparam value_t The underlying type
 * \ingroup search_configuration
 */
struct strata : detail::strong_type<uint8_t, strata, detail::strong_type_skill::convert>
{
    using detail::strong_type<uint8_t, strata, detail::strong_type_skill::convert>::strong_type;
};

/*!\brief Configuration element to determine the search mode.
 * \ingroup search_configuration
 */
template <typename mode_t>
//!\cond
    requires std::Same<remove_cvref_t<mode_t>, detail::search_mode_all> ||
             std::Same<remove_cvref_t<mode_t>, detail::search_mode_all_best> ||
             std::Same<remove_cvref_t<mode_t>, detail::search_mode_best> ||
             std::Same<remove_cvref_t<mode_t>, strata>
//!\endcond
class mode : public pipeable_config_element
{
public:

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::search_config_id id{detail::search_config_id::mode};

    //!\publicsecton
    /*!\name Constructor, destructor and assignment
     * \brief Defaulted all standard constructor.
     * \{
     */
    constexpr mode()                         noexcept = default;
    constexpr mode(mode const &)             noexcept = default;
    constexpr mode(mode &&)                  noexcept = default;
    constexpr mode & operator=(mode const &) noexcept = default;
    constexpr mode & operator=(mode &&)      noexcept = default;
    ~mode()                                  noexcept = default;

    /*!\brief Constructs an object from the given mode.
     * \tparam mode_t The type of the search mode.
     * \param  model  The mode to be used.
     */
    constexpr mode(mode_t mode) noexcept : value{std::move(mode)}
    {}
    //!}

    //!\brief The stored value.
    mode_t value{};
};

/*!\name Type deduction guides
 * \relates seqan3::search_cfg::mode
 * \{
 */
//!\brief Default type deduces to best mode.
mode() -> mode<detail::search_mode_best>;

//!\brief Deduces search mode type from constructor argument.
template <typename mode_t>
mode(mode_t) -> mode<remove_cvref_t<mode_t>>;
//!\}
} // namespace seqan3::search_cfg

//!\}
