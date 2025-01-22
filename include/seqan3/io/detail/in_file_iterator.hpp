// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides the seqan3::detail::in_file_iterator class template.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <cassert>
#include <ranges>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

/*!\brief Input iterator necessary for providing a range-like interface in input file.
 * \tparam file_type The data structure on which the iterator operates, e.g. `std::vector<int>`.
 * \implements std::input_Iterator
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
 * This iterator may be compared against std::default_sentinel_t, this check delegates to
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

    //!\brief The value type.
    using value_type = typename file_type::value_type;
    //!\brief The reference type.
    using reference = typename file_type::reference;
    //!\brief The const reference type.
    using const_reference = typename file_type::reference;
    //!\brief The size type.
    using size_type = typename file_type::size_type;
    //!\brief The difference type. A signed integer type, usually std::ptrdiff_t.
    using difference_type = typename file_type::difference_type;
    //!\brief The pointer type.
    using pointer = typename file_type::value_type *;
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
    constexpr in_file_iterator(in_file_iterator &&) = default;
    //!\brief Move assignment.
    constexpr in_file_iterator & operator=(in_file_iterator &&) = default;
    //!\brief Use default deconstructor.
    ~in_file_iterator() = default;

    //!\brief Construct with reference to host.
    constexpr in_file_iterator(file_type & _host) noexcept : host{&_host}
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
        return host->record_buffer;
    }

    //!\brief Dereference returns the currently buffered record.
    reference operator*() const noexcept
    {
        assert(host != nullptr);
        return host->record_buffer;
    }
    //!\}

    /*!\name Comparison operators
     * \brief Only (in-)equality comparison of iterator with end() is supported.
     * \{
     */

    //!\brief Checks whether `*this` is equal to the sentinel.
    constexpr bool operator==(std::default_sentinel_t const &) const noexcept
    {
        assert(host != nullptr);
        return host->at_end;
    }

    //!\brief Checks whether `*this` is not equal to the sentinel.
    constexpr bool operator!=(std::default_sentinel_t const &) const noexcept
    {
        assert(host != nullptr);
        return !host->at_end;
    }

    //!\brief Checks whether `it` is equal to the sentinel.
    constexpr friend bool operator==(std::default_sentinel_t const &, in_file_iterator const & it) noexcept
    {
        return (it == std::default_sentinel);
    }

    //!\brief Checks whether `it` is not equal to the sentinel.
    constexpr friend bool operator!=(std::default_sentinel_t const &, in_file_iterator const & it) noexcept
    {
        return (it != std::default_sentinel);
    }
    //!\}

    /*!\name File position functionality
     * \brief Low level API. Enables seeking to and returning specific file positions from the iterator.
     * \{
     */

    //!\brief Returns the current position in the file via `std::streampos`.
    std::streampos file_position() const
    {
        assert(host != nullptr);
        return host->position_buffer;
    }

    //!\brief Low level API. Sets the current position of the iterator to the given position, throws if at end of file.
    in_file_iterator & seek_to(std::streampos const & pos)
    {
        assert(host != nullptr);
        host->secondary_stream->seekg(pos);
        if (host->secondary_stream->fail())
        {
            throw std::runtime_error{"Seeking to file position failed!"};
        }
        host->at_end = false; // iterator will not be at end if seeking to a specific record
        host->read_next_record();
        return *this;
    }
    //!\}

private:
    //!\brief Pointer to file host.
    file_type * host{};
};

} // namespace seqan3::detail
