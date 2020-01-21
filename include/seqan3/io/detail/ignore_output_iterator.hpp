// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \brief Provides seqan3::detail::ignore_output_iterator for writing to `null` stream.
 */

#pragma once

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

/*!\brief An output iterator that emulates writing to a `null`-stream in order to dispose the output.
 * \ingroup io
 *
 * \details
 *
 * This helper iterator can be used to explicitly dispose output streams, i.e. the output data is transmitted nowhere.
 * A typical use case is when extracted bytes from an input stream should be ignored, as they contain only structural
 * information, e.g. a newline character.
 */
class ignore_output_iterator
{
public:

    /*!\name Member types
     * \{
     * \brief Associated types are void for output iterators, see also
     * [output iterator concept](https://en.cppreference.com/w/cpp/concept/output_iterator).
     */

    //!\brief The value type (void).
    using value_type        = void;
    //!\brief The reference type (void).
    using reference         = void;
    //!\brief The pointer type (void).
    using pointer           = void;
    //!\brief A signed integer type, usually std::ptrdiff_t.
    using difference_type   = std::ptrdiff_t;
    //!\brief The iterator category type.
    using iterator_category = std::output_iterator_tag;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    ignore_output_iterator() = default;                                            //!< Defaulted
    ignore_output_iterator(ignore_output_iterator const &) = default;              //!< Defaulted
    ignore_output_iterator(ignore_output_iterator &&) = default;                   //!< Defaulted
    ignore_output_iterator & operator= (ignore_output_iterator const &) = default; //!< Defaulted
    ignore_output_iterator & operator= (ignore_output_iterator &&) = default;      //!< Defaulted
    ~ignore_output_iterator() = default;                                           //!< Defaulted
    //!\}

    /*!\name Member functions
     * \brief Each function performs no operation. In fact writing to the seqan3::detail::ignore_output_iterator,
     *        is subject to removal by compiler optimizations.
     * \{
     */
    //!\brief Emulates writing the passed value to the `null`-stream.
    template <typename type>
    constexpr ignore_output_iterator & operator= (type const /*v*/) noexcept
    {
        return *this;
    }

    //!\brief This operator performs no function in output iterators.
    constexpr ignore_output_iterator & operator* () noexcept
    {
        return *this;
    }

    //!\brief This operator performs no function in output iterators.
    constexpr ignore_output_iterator & operator++ () noexcept
    {
        return *this;
    }

    //!\brief This operator performs no function in output iterators.
    constexpr ignore_output_iterator & operator++ (int) noexcept
    {
        return *this;
    }
    //!\}
};

} // namespace seqan3::detail
