// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::view::istreambuf.
 */

#pragma once

#include <iosfwd>

#include <seqan3/io/stream/concept.hpp>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

// ============================================================================
//  fast_istreambuf_iterator
// ============================================================================

//!\brief Functionally the same as std::basic_streambuf<char_t, traits_t_>, but exposes protected members as public.
template <typename char_t, typename traits_t = std::char_traits<char_t>>
struct stream_buffer_exposer : public std::basic_streambuf<char_t, traits_t>
{
    //!\brief The actual stream type.
    using base_t = std::basic_streambuf<char_t, traits_t>;

    // Expose protected members:
    using base_t::eback;
    using base_t::gptr;
    using base_t::egptr;
    using base_t::gbump;
    using base_t::underflow;

    using base_t::pbase;
    using base_t::pptr;
    using base_t::epptr;
    using base_t::pbump;
    using base_t::overflow;
};

/*!\brief Functionally the same as std::istreambuf_iterator, but faster.
 * \tparam char_t       The stream's character type.
 * \tparam traits_t_    The stream's traits type.
 *
 * \details
 *
 * Performs less virtual function calls than std::istreambuf_iterator.
 *
 * \todo Make this move-only after input iterators are allowed to be move-only.
 *
 */
template <typename char_t, typename traits_t = std::char_traits<char_t>>
class fast_istreambuf_iterator
{
private:
    //!\brief Down-cast pointer to the stream-buffer.
    stream_buffer_exposer<char_t, traits_t> * stream_buf = nullptr;

public:
    /*!\name Associated types
     * \{
     */
    using difference_type   = ptrdiff_t;               //!< Defaults to ptrdiff_t.
    using value_type        = char_t;                  //!< The char type of the stream.
    using reference         = char_t;                  //!< The char type of the stream.
    using pointer           = void;                    //!< Has no pointer type.
    using iterator_category = std::input_iterator_tag; //!< Pure input iterator.
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    fast_istreambuf_iterator()                                              noexcept = default; //!< Defaulted.
    fast_istreambuf_iterator(fast_istreambuf_iterator const &)              noexcept = default; //!< Defaulted.
    fast_istreambuf_iterator(fast_istreambuf_iterator &&)                   noexcept = default; //!< Defaulted.
    fast_istreambuf_iterator & operator=(fast_istreambuf_iterator const &)  noexcept = default; //!< Defaulted.
    fast_istreambuf_iterator & operator=(fast_istreambuf_iterator &&)       noexcept = default; //!< Defaulted.
    ~fast_istreambuf_iterator()                                             noexcept = default; //!< Defaulted.

    //!\brief Construct from a stream buffer.
    explicit fast_istreambuf_iterator(std::basic_streambuf<char_t, traits_t> & ibuf) :
        stream_buf{reinterpret_cast<stream_buffer_exposer<char_t, traits_t> *>(&ibuf)}
    {
        assert(stream_buf != nullptr);
        stream_buf->underflow(); // ensure the stream buffer has content on construction
    }
    //!\}

    /*!\name Arithmetic operators
     * \{
     */
    //!\brief Advance by one and rebuffer if necessary (vtable lookup iff rebuffering).
    fast_istreambuf_iterator & operator++()
    {
        assert(stream_buf != nullptr);
        stream_buf->snextc(); // move to then right, then underflow()
        return *this;
    }

    //!\overload
    void operator++(int)
    {
        ++(*this);
    }
    //!\}

    //!\brief Read current value from buffer (no vtable lookup, safe if not at end).
    reference operator*() const
    {
        assert(stream_buf != nullptr);
        return *stream_buf->gptr();
    }

    /*!\name Comparison operators
     * \brief We define comparison only against the sentinel.
     * \{
     */
    //!\brief True if the read buffer is not empty; involves no vtable lookup.
    friend bool operator==(fast_istreambuf_iterator const & lhs, std::ranges::default_sentinel_t const &) noexcept
    {
        assert(lhs.stream_buf != nullptr);
        // compare size of remaining buffer; since ++ always resizes if possible, safe to compare pointers here
        return (lhs.stream_buf->gptr() == lhs.stream_buf->egptr());
    }

    //!\brief True if the read buffer is empty; involves no vtable lookup.
    friend bool operator!=(fast_istreambuf_iterator const & lhs, std::ranges::default_sentinel_t const &) noexcept
    {
        return !(lhs == std::ranges::default_sentinel);
    }

    //!\brief True if the read buffer is not empty; involves no vtable lookup.
    friend bool operator==(std::ranges::default_sentinel_t const &, fast_istreambuf_iterator const & rhs) noexcept
    {
        return rhs == std::ranges::default_sentinel;
    }

    //!\brief True if the read buffer is empty; involves no vtable lookup.
    friend bool operator!=(std::ranges::default_sentinel_t const &, fast_istreambuf_iterator const & rhs) noexcept
    {
        return !(rhs == std::ranges::default_sentinel);
    }
    //!\}
};

// ============================================================================
//  istreambuf_fn (adaptor definition)
// ============================================================================

//!\brief View adaptor/factory definition for view::istream.
struct istreambuf_fn
{
    /*!\brief Return the view object.
     * \param[in,out] s Reference to the stream buffer.
     * \tparam stream_char_t Character type of the stream device.
     * \tparam stream_traits_t Traits type of the stream device.
     * \returns A std::ranges::subrange over a detail::fast_istreambuf_iterator and std::ranges::default_sentinel_t.
     */
    template <typename stream_char_t, typename stream_traits_t>
    constexpr auto operator()(std::basic_streambuf<stream_char_t, stream_traits_t> & s) const
    {
        return std::ranges::subrange<detail::fast_istreambuf_iterator<stream_char_t, stream_traits_t>,
                                     std::ranges::default_sentinel_t>
        {
            detail::fast_istreambuf_iterator<stream_char_t, stream_traits_t>{s},
            std::ranges::default_sentinel_t{}
        };
    }

    /*!\brief Return the view object.
     * \tparam stream_t Type of the stream, must model seqan3::IStream2.
     * \param[in,out] s Reference to a stream object.
     * \returns A std::ranges::subrange over a detail::fast_istreambuf_iterator and std::ranges::default_sentinel_t.
     */
    template <IStream2 stream_t>
    constexpr auto operator()(stream_t & s) const
    {
        return this->operator()(*s.rdbuf());
    }
};

} // namespace seqan3::detail

// ============================================================================
//  view::istreambuf (adaptor instance definition)
// ============================================================================

namespace seqan3::view
{

/*!\name General purpose views
 * \{
 */

/*!\brief                A view factory that returns a view over the stream buffer of an input stream.
 * \tparam istreambuf_t  The type of the stream(buffer); must be std::basic_streambuf or model seqan3::IStream2.
 * \param[in] istreambuf The stream buffer or an input stream of whome the buffer is retrieved.
 * \returns
 * \ingroup view
 *
 * \details
 *
 * **Header**
 * ```cpp
 *      #include <seqan3/range/view/istreambuf.hpp>
 * ```
 *
 * ### View properties
 *
 * This is a source-only view adaptor, also known as a range factory; you cannot pipe anything into it.
 *
 * | range concepts and reference_t  | `rrng_t` (returned range type)   |
 * |---------------------------------|:--------------------------------:|
 * | std::ranges::InputRange         | *guaranteed*                     |
 * | std::ranges::ForwardRange       |                                  |
 * | std::ranges::BidirectionalRange |                                  |
 * | std::ranges::RandomAccessRange  |                                  |
 * | std::ranges::ContiguousRange    |                                  |
 * |                                 |                                  |
 * | std::ranges::ViewableRange      | *guaranteed*                     |
 * | std::ranges::View               | *guaranteed*                     |
 * | std::ranges::SizedRange         |                                  |
 * | std::ranges::CommonRange        |                                  |
 * | std::ranges::OutputRange        |                                  |
 * | seqan3::ConstIterableRange      | *guaranteed*                     |
 * |                                 |                                  |
 * | seqan3::reference_t             | `istream_t::char_type`           |
 *
 * See the \link view view submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * This adaptor is different from std::ranges::istream_range in that it operates directly on the buffer.
 * It further uses a custom streambuf_iterator (not std::istreambuf_iterator) that performs less virtual
 * function calls.
 *
 * \hideinitializer
 */
inline constexpr auto istreambuf = detail::istreambuf_fn{};
//!\}

} // namespace seqan3::view
