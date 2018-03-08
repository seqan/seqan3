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
 * \brief Provides the seqan3::detail::in_file_iterator class.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <cassert>

#include <seqan3/core/platform.hpp>
#include <seqan3/range/concept.hpp>

namespace seqan3::detail
{

template <typename file_type>
class in_file_sentinel
{};


/*!\brief Implementation of a random access iterator on an input container pointer.
 * \tparam file_type The data structure on which the iterator operates, e.g. `std::vector<int>`.
 *
 */
template <typename file_type>
class in_file_iterator
{
public:
    //!\brief The value_type is the recored_type.
    using value_type = typename file_type::value_type;
    //!\brief The reference_type is an rvalue reference, because the record is always moved out..
    using reference = typename file_type::reference;
    //!\brief The const_reference type.
    using const_reference = typename file_type::const_reference;
    //!\brief An unsigned integer type, usually std::size_t.
    using size_type = typename file_type::size_type;
    //!\brief A signed integer type, usually std::ptrdiff_t.
    using difference_type = typename file_type::difference_type;

    //!\brief Tag this class as an input access iterator.
    using iterator_category = std::input_iterator_tag;

    /*!\name Constructors/Destructors
     * \{
    */
    //!\brief Default constructor.
    constexpr in_file_iterator() = default;
    //!\brief Construct by host, default position pointer with 0.
    constexpr in_file_iterator(file_type & _host) noexcept :
        host{&_host}
    {}
    //!\brief Copy constructor.
    constexpr in_file_iterator(in_file_iterator const &) = default;
    //!\brief Copy construction via assignment.
    constexpr in_file_iterator & operator=(in_file_iterator const &) = default;
    //!\brief Move constructor.
    constexpr in_file_iterator (in_file_iterator &&) = default;
    //!\brief Move assignment.
    constexpr in_file_iterator & operator=(in_file_iterator &&) = default;
    //!\brief Use default deconstructor.
    ~in_file_iterator() = default;

    //!\brief Constructor for const version from non-const version.
    template <typename file_type2>
    //!\cond
        requires std::is_const_v<file_type> && !std::is_const_v<file_type2> &&
                 std::is_same_v<std::remove_const_t<file_type>, file_type2>
    //!\endcond
    constexpr in_file_iterator(in_file_iterator<file_type2> const & rhs) noexcept :
        host{rhs.host}
    {}
    //!\}

    template <typename = std::enable_if_t<!std::is_const_v<file_type>>>
    in_file_iterator & operator++()
    {
        assert(host != nullptr);
        host->buffer_next_record();
        return *this;
    }

    //!\brief Post-increment, return previous iterator state. TODO warning previous is "invalid"
    template <typename = std::enable_if_t<!std::is_const_v<file_type>>>
    in_file_iterator operator++(int)
    {
        assert(host != nullptr);
        in_file_iterator cpy{*this};
        host->buffer_next_record();
        return cpy;
     }

     template <typename = std::enable_if_t<!std::is_const_v<file_type>>>
     reference operator*() const noexcept
     {
         return std::move(host->back());
     }

    /*!\name Comparison operators
     * \brief Compares only the absolute position of two iterators.
     * \{
     */
    constexpr bool operator==(in_file_sentinel<file_type> const &) const noexcept
    {
        assert(host != nullptr);
        return host->stream.eof();
    }

    constexpr bool operator!=(in_file_sentinel<file_type> const &) const noexcept
    {
        assert(host != nullptr);
        return !host->stream.eof();
    }

    constexpr friend bool operator==(in_file_sentinel<file_type> const &,
                                     in_file_iterator const & it) noexcept
    {
        assert(host != nullptr);
        return (it == in_file_sentinel<file_type>{});
    }

    constexpr friend bool operator!=(in_file_sentinel<file_type> const &,
                                     in_file_iterator const & it) noexcept
    {
        assert(host != nullptr);
        return (it != in_file_sentinel<file_type>{});
    }
    //!\}
private:
    file_type * host = nullptr;

    //!\brief The non-const version needs to give the const version access.
    template <typename file_type2>
    //!\cond
        requires std::is_const_v<file_type2> && !std::is_const_v<file_type> &&
                 std::is_same_v<std::remove_const_t<file_type2>, file_type>
    //!\endcond
    friend class in_file_iterator;
};

} // namespace seqan3::detail
