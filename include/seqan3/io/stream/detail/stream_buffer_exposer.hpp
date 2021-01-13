// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::stream_buffer_exposer.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <iosfwd>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

/*!\brief Functionally the same as std::basic_streambuf<char_t, traits_t_>, but exposes protected members as public.
 * \ingroup stream
 * \tparam char_t   The stream's character type.
 * \tparam traits_t The stream's traits type.
 *
 * \details
 *
 * This wrapper adds no functionality to std::basic_streambuf and is only used to expose protected members to
 * access the get and put area of the std::basic_streambuf.
 */
template <typename char_t, typename traits_t = std::char_traits<char_t>>
struct stream_buffer_exposer : public std::basic_streambuf<char_t, traits_t>
{
    //!\brief The actual stream type.
    using base_t = std::basic_streambuf<char_t, traits_t>;

    //!\cond
    // Expose protected members:
    using base_t::eback;
    using base_t::gptr;
    using base_t::egptr;
    using base_t::setg;
    using base_t::gbump;
    using base_t::underflow;

    using base_t::pbase;
    using base_t::pptr;
    using base_t::epptr;
    using base_t::pbump;
    using base_t::overflow;
    //!\endcond
};
} // namespace seqan3::detail
