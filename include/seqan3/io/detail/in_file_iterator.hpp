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
 * \brief Provides the seqan3::detail::in_file_iterator class template.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <cassert>

#include <range/v3/range_fwd.hpp>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

/*!\brief Input iterator necessary for providing a range-like interface in input file.
 * \tparam file_type The data structure on which the iterator operates, e.g. `std::vector<int>`.
 * \implements seqan3::input_iterator_concept
 * \ingroup io
 *
 * This iterator is a single-pass input iterator for input files. All member types are resolved
 * via `file_type`'s member types, dereference is implemented via file's `front()` member
 * function and increment calls the `buffer_next_record()` member of file.
 *
 * Note that since this is a single-pass input iterator, post-increment returns void because
 * previous iterators are always invalid (all iterators point to the current position in single-pass
 * ranges).
 *
 * This iterator may be compared against ranges::default_sentinel, this check delegates to
 * calling the `eof()` member function on the file's stream.
 */
template <typename file_type>
class in_file_iterator
{
    static_assert(!std::is_const_v<file_type>,
                  "You cannot iterate over const files, because the iterator changes the file.");
public:
    /*!\name Member types
     * \brief The associated types are derived from the `file_type`.
     * \{
     */
    using value_type        = typename file_type::value_type;
    using reference         = typename file_type::reference;
    using const_reference   = typename file_type::reference;
    using size_type         = typename file_type::size_type;
    using difference_type   = typename file_type::difference_type;
    //!\brief Tag this class as an input iterator.
    using iterator_category = std::input_iterator_tag;
    //!\}

    /*!\name Constructors, destructor and assignment.
     * \{
     */
    //!\brief Default constructor.
    constexpr in_file_iterator() = default;
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

    //!\brief Construct with reference to host.
    constexpr in_file_iterator(file_type & _host) noexcept :
        host{&_host}
    {}
    //!\}

    /*!\name Iterator operations
     * \{
     */
    //!\brief Move to the next record in the file and return a reference to it.
    in_file_iterator & operator++()
    {
        assert(host != nullptr);
        host->read_next_record();
        return *this;
    }

    //!\brief Post-increment is the same as pre-increment, but returns void.
    void operator++(int)
    {
        assert(host != nullptr);
        ++(*this);
    }

    //!\brief Dereference returns the currently buffered record.
    reference operator*() noexcept
    {
        assert(host != nullptr);
        return host->front();
    }

    //!\brief Dereference returns the currently buffered record.
    reference operator*() const noexcept
    {
        assert(host != nullptr);
        return host->front();
    }
    //!\}

    /*!\name Comparison operators
     * \brief Only (in-)equality comparison of iterator with end() is supported.
     * \{
     */
    constexpr bool operator==(ranges::default_sentinel const &) const noexcept
    {
        assert(host != nullptr);
        return host->stream.eof();
    }

    constexpr bool operator!=(ranges::default_sentinel const &) const noexcept
    {
        assert(host != nullptr);
        return !host->stream.eof();
    }

    constexpr friend bool operator==(ranges::default_sentinel const &,
                                     in_file_iterator const & it) noexcept
    {
        return (it == ranges::default_sentinel{});
    }

    constexpr friend bool operator!=(ranges::default_sentinel const &,
                                     in_file_iterator const & it) noexcept
    {
        return (it != ranges::default_sentinel{});
    }
    //!\}

private:
    //!\brief Pointer to file host.
    file_type * host{};
};

} // namespace seqan3::detail
