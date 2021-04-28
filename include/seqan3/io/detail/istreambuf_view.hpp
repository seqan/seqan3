// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::detail::istreambuf.
 */

#pragma once

#include <iosfwd>

#include <seqan3/std/iterator>
#include <seqan3/std/ranges>

#include <seqan3/io/stream/concept.hpp>
#include <seqan3/io/stream/detail/fast_istreambuf_iterator.hpp>

namespace seqan3::detail
{

// ============================================================================
//  istreambuf_fn (adaptor definition)
// ============================================================================

//!\brief View adaptor/factory definition for views::istream.
struct istreambuf_fn
{
    /*!\brief Return the view object.
     * \param[in,out] s Reference to the stream buffer.
     * \tparam stream_char_t builtin_characteracter type of the stream device.
     * \tparam stream_traits_t Traits type of the stream device.
     * \returns A std::ranges::subrange over a detail::fast_istreambuf_iterator and std::default_sentinel_t.
     */
    template <typename stream_char_t, typename stream_traits_t>
    constexpr auto operator()(std::basic_streambuf<stream_char_t, stream_traits_t> & s) const
    {
        return std::ranges::subrange<detail::fast_istreambuf_iterator<stream_char_t, stream_traits_t>,
                                     std::default_sentinel_t>
        {
            detail::fast_istreambuf_iterator<stream_char_t, stream_traits_t>{s},
            std::default_sentinel_t{}
        };
    }

    /*!\brief Return the view object.
     * \tparam stream_t Type of the stream, must model seqan3::input_stream.
     * \param[in,out] s Reference to a stream object.
     * \returns A std::ranges::subrange over a detail::fast_istreambuf_iterator and std::default_sentinel_t.
     */
    template <input_stream stream_t>
    constexpr auto operator()(stream_t & s) const
    {
        return this->operator()(*s.rdbuf());
    }
};

} // namespace seqan3::detail

// ============================================================================
//  detail::istreambuf (adaptor instance definition)
// ============================================================================

namespace seqan3::detail
{

/*!\name General purpose views
 * \{
 */

//!\copydoc seqan3::views::istreambuf
inline constexpr auto istreambuf = detail::istreambuf_fn{};
//!\}

} // namespace seqan3::detail

#ifdef SEQAN3_DEPRECATED_310
namespace seqan3::views
{
/*!\brief                A view factory that returns a view over the stream buffer of an input stream.
 * \tparam istreambuf_t  The type of the stream(buffer); must be std::basic_streambuf or model seqan3::input_stream.
 * \param[in] istreambuf The stream buffer or an input stream of whome the buffer is retrieved.
 * \returns
 * \ingroup views
 *
 * \details
 *
 * \header_file{seqan3/io/detail/istreambuf_view.hpp}
 *
 * ### View properties
 *
 * This is a source-only view adaptor, also known as a range factory; you cannot pipe anything into it.
 *
 * | Concepts and traits              | `rrng_t` (returned range type)   |
 * |----------------------------------|:--------------------------------:|
 * | std::ranges::input_range         | *guaranteed*                     |
 * | std::ranges::forward_range       |                                  |
 * | std::ranges::bidirectional_range |                                  |
 * | std::ranges::random_access_range |                                  |
 * | std::ranges::contiguous_range    |                                  |
 * |                                  |                                  |
 * | std::ranges::viewable_range      | *guaranteed*                     |
 * | std::ranges::view                | *guaranteed*                     |
 * | std::ranges::sized_range         |                                  |
 * | std::ranges::common_range        |                                  |
 * | std::ranges::output_range        |                                  |
 * | seqan3::const_iterable_range     | *guaranteed*                     |
 * |                                  |                                  |
 * | std::ranges::range_reference_t   | `istream_t::char_type`           |
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * This adaptor is different from std::ranges::istream_range in that it operates directly on the buffer.
 * It further uses a custom streambuf_iterator (not std::istreambuf_iterator) that performs less virtual
 * function calls.
 *
 * \hideinitializer
 *
 * \deprecated Use std::ranges::istream_view(std::basic_istream) instead (you need to pass a stream instead of a
 *             std::basic_streambuf)
 */
SEQAN3_DEPRECATED_310 inline constexpr auto istreambuf = detail::istreambuf;
} // namespace seqan3::views
#endif // SEQAN3_DEPRECATED_310
