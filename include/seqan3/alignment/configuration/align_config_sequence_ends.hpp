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
 * \brief Provides the free ends gap configuration.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

#pragma once

#include <iostream>

#include <seqan3/alignment/configuration/utility.hpp>
#include <seqan3/core/add_enum_bitwise_operators.hpp>
#include <seqan3/core/algorithm/all.hpp>
#include <seqan3/core/metafunction/template_inspection.hpp>

namespace seqan3
{

/*!\brief Enum class for all supported sequence ends specifications.
 * \ingroup configuration
 * \details The members specify where continuous gaps in the beginning or end of a sequence are not penalized
 * in the alignment.
 */
enum struct free_ends_at : uint8_t
{
    //!\brief No free gaps at the sequence ends. Each gap is scored according to seqan3::align_cfg::gap.
    none       = 0b0000,
    //!\brief Continuous gaps in the beginning of the first sequence are not scored.
    seq1_front = 0b0001,
    //!\brief Continuous gaps at the end of the first sequence are not scored.
    seq1_back  = 0b0010,
    //!\brief Continuous gaps in the beginning of the second sequence are not scored.
    seq2_front = 0b0100,
    //!\brief Continuous gaps at the end of the second sequence are not scored.
    seq2_back  = 0b1000,
    //!\brief Continuous gaps in the beginning and end of the first sequence are not scored.
    seq1       = seq1_front | seq1_back,
    //!\brief Continuous gaps in the beginning and end of the second sequence are not scored.
    seq2       = seq2_front | seq2_back,
    //!\brief Continuous gaps in the beginning and end of both sequences are not scored.
    all        = seq1 | seq2
};

//!\cond
template <>
constexpr bool add_enum_bitwise_operators<seqan3::free_ends_at> = true;
//!\endcond

} // namespace seqan3

namespace seqan3::detail
{
/*!\brief A configuration element for gaps at the sequence ends.
 * \ingroup configuration
 * \tparam val An enum value that specifies which sequence ends allow gaps without penalty.
 */
template <free_ends_at val = free_ends_at::none>
struct align_config_sequence_ends
{
    //!\brief Holds the actual setting.
    static constexpr free_ends_at value = val;
};

/*!\brief A deferred configuration element for gaps at the sequence ends.
 * \ingroup configuration
 */
struct align_config_sequence_ends_deferred : public detail::deferred_config_element_base<align_config_sequence_ends_deferred>
{
    //!\brief Holds the actual setting.
    free_ends_at value{};

    /*!
     * \brief Adds to the configuration a configuration element for free gaps at the sequence ends.
     * \tparam fn_t Type of the function to be invoked.
     * \tparam configuration_t The type of the configuration to be extended.
     * \param fn Callable to be invoked for the configuration.
     * \param cfg The configuration to be extended.
     * \returns A new configuration containing the sequence ends configuration element.
     */
    template <typename fn_t, typename configuration_t>
    constexpr auto invoke(fn_t && fn, configuration_t && cfg) const
    //!\cond
        requires detail::is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
    //!\endcond
    {
        switch (static_cast<uint8_t>(value))
        {
            case 0b0000: return fn(std::forward<configuration_t>(cfg).replace_with(*this,
                                   align_config_sequence_ends<static_cast<free_ends_at>(0b0000)>{}));

            case 0b0001: return fn(std::forward<configuration_t>(cfg).replace_with(*this,
                                   align_config_sequence_ends<static_cast<free_ends_at>(0b0001)>{}));

            case 0b0010: return fn(std::forward<configuration_t>(cfg).replace_with(*this,
                                   align_config_sequence_ends<static_cast<free_ends_at>(0b0010)>{}));

            case 0b0011: return fn(std::forward<configuration_t>(cfg).replace_with(*this,
                                   align_config_sequence_ends<static_cast<free_ends_at>(0b0011)>{}));

            case 0b0100: return fn(std::forward<configuration_t>(cfg).replace_with(*this,
                                   align_config_sequence_ends<static_cast<free_ends_at>(0b0100)>{}));

            case 0b0101: return fn(std::forward<configuration_t>(cfg).replace_with(*this,
                                   align_config_sequence_ends<static_cast<free_ends_at>(0b0101)>{}));

            case 0b0110: return fn(std::forward<configuration_t>(cfg).replace_with(*this,
                                   align_config_sequence_ends<static_cast<free_ends_at>(0b0110)>{}));

            case 0b0111: return fn(std::forward<configuration_t>(cfg).replace_with(*this,
                                   align_config_sequence_ends<static_cast<free_ends_at>(0b0111)>{}));

            case 0b1000: return fn(std::forward<configuration_t>(cfg).replace_with(*this,
                                   align_config_sequence_ends<static_cast<free_ends_at>(0b1000)>{}));

            case 0b1001: return fn(std::forward<configuration_t>(cfg).replace_with(*this,
                                   align_config_sequence_ends<static_cast<free_ends_at>(0b1001)>{}));

            case 0b1010: return fn(std::forward<configuration_t>(cfg).replace_with(*this,
                                   align_config_sequence_ends<static_cast<free_ends_at>(0b1010)>{}));

            case 0b1011: return fn(std::forward<configuration_t>(cfg).replace_with(*this,
                                   align_config_sequence_ends<static_cast<free_ends_at>(0b1011)>{}));

            case 0b1100: return fn(std::forward<configuration_t>(cfg).replace_with(*this,
                                   align_config_sequence_ends<static_cast<free_ends_at>(0b1100)>{}));

            case 0b1101: return fn(std::forward<configuration_t>(cfg).replace_with(*this,
                                   align_config_sequence_ends<static_cast<free_ends_at>(0b1101)>{}));

            case 0b1110: return fn(std::forward<configuration_t>(cfg).replace_with(*this,
                                   align_config_sequence_ends<static_cast<free_ends_at>(0b1110)>{}));

            case 0b1111: return fn(std::forward<configuration_t>(cfg).replace_with(*this,
                                   align_config_sequence_ends<static_cast<free_ends_at>(0b1111)>{}));

            default:     throw std::invalid_argument("Enum value out of bounds for seqan3::free_ends_at.");
        }
    }
};

/*!\brief The sequence ends adaptor enabling pipe notation.
 * \ingroup configuration
 * \tparam val An enum value that specifies which sequence ends allow gaps without penalty.
 */
template <free_ends_at val>
struct align_config_sequence_ends_adaptor : public configuration_fn_base<align_config_sequence_ends_adaptor<val>>
{
    /*!\brief Adds to the configuration a configuration element for free gaps at the sequence ends.
     * \tparam configuration_type The type of the configuration to be extended.
     * \param[in] cfg The configuration to be extended.
     * \param[in] _val The value to be set.
     * \returns A new configuration containing the sequence ends configuration element.
     */
    template <typename configuration_type>
    //!\cond
        requires is_algorithm_configuration_v<remove_cvref_t<configuration_type>>
    //!\endcond
    constexpr auto invoke(configuration_type && cfg, free_ends_at const _val) const
    {
        static_assert(is_valid_alignment_configuration_v<align_cfg::id::sequence_ends, remove_cvref_t<configuration_type>>,
                      SEQAN3_INVALID_CONFIG(align_cfg::id::sequence_ends));

        align_config_sequence_ends_deferred tmp;
        tmp.value = _val;
        return std::forward<configuration_type>(cfg).push_front(std::move(tmp));
    }

    /*!\brief Adds to the configuration a configuration element for free gaps at the sequence ends.
     * \tparam configuration_type The type of the configuration to be extended.
     * \param[in] cfg The configuration to be extended.
     * \returns A new configuration containing the sequence ends configuration element.
     * \details This function uses the template value.
     */
    template <typename configuration_type>
    //!\cond
        requires is_algorithm_configuration_v<remove_cvref_t<configuration_type>>
    //!\endcond
    constexpr auto invoke(configuration_type && cfg) const
    {
        static_assert(is_valid_alignment_configuration_v<align_cfg::id::sequence_ends, remove_cvref_t<configuration_type>>,
                      SEQAN3_INVALID_CONFIG(align_cfg::id::sequence_ends));

        return std::forward<configuration_type>(cfg).push_front(align_config_sequence_ends<val>{});
    }
};

//!\brief Helper template meta-function associated with seqan3::detail::align_config_sequence_ends.
//!\ingroup configuration
template <>
struct on_align_config<align_cfg::id::sequence_ends>
{
    /*!\brief Type alias used by meta::find_if
     * \tparam cfg The configuration element.
     */
    template <config_element_concept cfg>
    using invoke = std::conditional_t<is_value_specialisation_of_v<cfg, align_config_sequence_ends>,
                                      std::true_type,
                                      typename std::is_same<cfg, align_config_sequence_ends_deferred>::type>;
};

/*!\brief Mapping from the seqan3::detail::align_config_sequence_ends type to its corresponding seqan3::align_cfg::id.
 * \ingroup configuration
 * \tparam val An enum value that specifies which sequence ends allow gaps without penalty.
 */
template <free_ends_at val>
struct align_config_type_to_id<align_config_sequence_ends<val>>
{
    //!\brief The associated seqan3::align_cfg::id.
    static constexpr align_cfg::id value = align_cfg::id::sequence_ends;
};

/*!\brief Mapping from the seqan3::detail::align_config_sequence_ends type to its corresponding seqan3::align_cfg::id.
 * \ingroup configuration
 */
template <>
struct align_config_type_to_id<align_config_sequence_ends_deferred>
{
    //!\brief The associated seqan3::align_cfg::id.
    static constexpr align_cfg::id value = align_cfg::id::sequence_ends;
};

} // namespace seqan3::detail

namespace seqan3::align_cfg
{

/*!\brief A configuration adaptor for gaps at the sequence ends.
 * \ingroup configuration
 * \tparam val An enum value that specifies which sequence ends allow gaps without penalty.
 *
 * \details This configuration allows to specify, whether continuous gaps in the front or end of a sequence
 * are penalized in the alignment.
 */
template <free_ends_at val = free_ends_at::none>
inline constexpr detail::align_config_sequence_ends_adaptor<val> sequence_ends;

} // namespace seqan3::align_cfg
