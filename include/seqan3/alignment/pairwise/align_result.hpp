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
 * \brief Provides seqan3::align_result.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <functional>
#include <optional>
#include <tuple>

#include <seqan3/core/metafunction/template_inspection.hpp>

namespace seqan3
{

/*!\brief Keys for different alignment results.
 * \ingroup pairwise
 */
enum struct align_result_key : uint8_t
{
    id,     // cannot be set by the config
    score,  // report the score
    end,    // report the end position and score
    begin,  // report the begin and end position and score
    trace   // report the full trace and begin and end and score
};

/*!\brief Stores the alignment results and offers a tuple-like interface.
 * \extends std::tuple
 * \ingroup pairwise
 * \tparam output_type_list_t seqan3::type_list over the contained results.
 */
template <typename output_type_list_t>
struct align_result //!\cond
    : public detail::transfer_template_args_onto_t<output_type_list_t, std::tuple>
    //!\endcond
{

    static_assert(meta::size<output_type_list_t>::value >= 2,
                  "Logic error: This class requires at least a type list with two template arguments "
                  "for the id and the alignment score");

    //!\brief The base tuple type.
    using base_type = detail::transfer_template_args_onto_t<output_type_list_t, std::tuple>;

    //!\brief Inheriting the constructors of the base type.
    using base_type::base_type;

    /*!\name Access functions
     * \brief Convenience helper functions to access elements of the alignment result type.
     * \{
     */

    //!\brief Returns the id.
    constexpr std::tuple_element_t<0, base_type> id() const noexcept
    {
        return std::get<0>(*this);
    }

    //!\brief Returns the score.
    constexpr std::tuple_element_t<1, base_type> score() const noexcept
    {
        return std::get<1>(*this);
    }

    /*!\brief Returns the end coordinate of the alignment if requested by the algorithm configuration,
     *        otherwise `std::ignore`.
     */
    constexpr decltype(auto) end_coordinate() const & noexcept
    {
        if constexpr (std::tuple_size_v<base_type> < 3)
            return std::ignore;
        else
            return std::get<2>(*this);
    }

    //!\overload
    constexpr decltype(auto) end_coordinate() const && noexcept
    {
        if constexpr (std::tuple_size_v<base_type> < 3)
            return std::ignore;
        else //TODO Fix outer move after gcc-7
            return std::move(std::get<2>(std::move(*this)));
    }

    /*!\brief Returns the begin coordinate of the alignment if requested by the algorithm configuration,
     *        otherwise `std::ignore`.
     */
    constexpr decltype(auto) begin_coordinate() const & noexcept
    {
        if constexpr (std::tuple_size_v<base_type> < 4)
            return std::ignore;
        else
            return std::get<3>(*this);
    }

    //!\overload
    constexpr decltype(auto) begin_coordinate() const && noexcept
    {
        if constexpr (std::tuple_size_v<base_type> < 4)
            return std::ignore;
        else //TODO Fix outer move after gcc-7 fix
            return std::move(std::get<3>(std::move(*this)));
    }

    /*!\brief Returns the traceback of the alignment if requested by the algorithm configuration,
     *        otherwise `std::ignore`.
     */
    constexpr decltype(auto) trace() const & noexcept
    {
        if constexpr (std::tuple_size_v<base_type> < 5)
            return std::ignore;
        else
            return std::get<4>(*this);
    }

    //!\overload
    constexpr decltype(auto) trace() const && noexcept
    {
        if constexpr (std::tuple_size_v<base_type> < 5)
            return std::ignore;
        else //TODO Fix outer move after gcc-7 fix
            return std::move(std::get<4>(std::move(*this)));
    }
     //!\}
};

/*!\name Tuple-like get interface
 * \ingroup pairwise
 * \relates seqan3::align_result
 * \{
 */

/*!\brief Returns the specified element using seqan3::align_result_key as specifier.
 * \tparam    e         The enum value of seqan3::align_result_key.
 * \param[in] align_res The seqan3::align_result to get the specified value for.
 * \returns The value associated with the enum value `e`.
 */
template <align_result_key e, typename output_type_list_t>
constexpr auto & get(align_result<output_type_list_t> & align_res) noexcept
{
    constexpr size_t index = static_cast<uint8_t>(e);

    static_assert(index < std::tuple_size_v<align_result<output_type_list_t>>,
                  "Element access error: index out of range.");
    return std::get<index>(align_res);
}

//!\overload
template <align_result_key e, typename output_type_list_t>
constexpr auto const & get(align_result<output_type_list_t> const & align_res) noexcept
{
    constexpr size_t index = static_cast<uint8_t>(e);

    static_assert(index < std::tuple_size_v<align_result<output_type_list_t>>,
                  "Element access error: index out of range.");
    return std::get<index>(align_res);
}

//!\overload
template <align_result_key e, typename output_type_list_t>
constexpr auto && get(align_result<output_type_list_t> && align_res) noexcept
{
    constexpr size_t index = static_cast<uint8_t>(e);

    static_assert(index < std::tuple_size_v<align_result<output_type_list_t>>,
                  "Element access error: index out of range.");
    return std::get<index>(std::move(align_res));
}

//!\overload
template <align_result_key e, typename output_type_list_t>
constexpr auto const && get(align_result<output_type_list_t> const && align_res) noexcept
{
    constexpr size_t index = static_cast<uint8_t>(e);

    static_assert(index < std::tuple_size_v<align_result<output_type_list_t>>,
                  "Element access error: index of range.");

    // TODO Remove superfluous outer move after fixed in gcc-7
    return std::move(std::get<index>(std::move(align_res)));
}
//!\}
} // namespace seqan3

namespace std
{

/*!\brief Overloads the tuple element type trait function for seqan3::align_result.
 * \ingroup pairwise
 */
template <size_t idx, typename output_type_list_t>
struct tuple_element<idx, seqan3::align_result<output_type_list_t>>
{
    //!\brief The element type at the given position.
    using type = std::tuple_element_t<idx, typename seqan3::align_result<output_type_list_t>::base_type>;
};

/*!\brief Overloads the tuple size type trait function for seqan3::align_result.
 * \ingroup pairwise
 */
template <typename output_type_list_t>
struct tuple_size<seqan3::align_result<output_type_list_t>>
{
    //!\brief The number of elements.
    static constexpr size_t value = std::tuple_size_v<typename seqan3::align_result<output_type_list_t>::base_type>;
};
} // namespace std
