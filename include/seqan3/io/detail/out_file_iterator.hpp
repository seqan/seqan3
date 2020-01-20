// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the seqan3::detail::out_file_iterator class template.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <cassert>

#include <range/v3/range_fwd.hpp>

#include <seqan3/core/platform.hpp>

#include <seqan3/std/ranges>

namespace seqan3::detail
{

/*!\brief Output iterator necessary for providing a range-like interface in output file.
 * \tparam file_type The data structure on which the iterator operates, e.g. `std::vector<int>`.
 * \implements std::output_Iterator
 * \ingroup io
 *
 * This iterator is a single-pass output iterator for output files. All associated types
 * are void, dereferencing this iterator returns a reference to self so that you can use
 * the assignment operator, both, on the iterator directly, and after dereferencing.
 * The assignment operator performs the actual file-writing, it delegates to the file's
 * `push_back()` member functions. The argument given to `push_back()` must meet the
 * specific file's constraints and is not checked by the iterator.
 *
 * The increment operations are no-ops (perform nothing) and comparisons against
 * std::ranges::default_sentinel always return false (there is no end in an output file).
 *
 * If any of these characteristics seem unusual to you, please refer to the [standard library's
 * documentation on output iterators](https://en.cppreference.com/w/cpp/concept/output_iterator).
 *
 * This class template differs from std::back_insert_iterator only in that it performs
 * no checks itself on the assigned values and that it allows comparisons against
 * std::ranges::default_sentinel.
 */
template <typename file_type>
class out_file_iterator
{
    static_assert(!std::is_const_v<file_type>,
                  "You cannot iterate over const files, because the iterator changes the file.");
public:
    /*!\name Member types
     * \brief The associated types are `void`, see the full description.
     * \{
     */

    //!\brief The value type (void).
    using value_type        = void;
    //!\brief The reference type (void).
    using reference         = void;
    //!\brief The const reference type (void).
    using const_reference   = void;
    //!\brief The size type (void).
    using size_type         = void;
    //!\brief A signed integer type, usually std::ptrdiff_t.
    using difference_type   = std::ptrdiff_t;
    //!\brief The pointer type.
    using pointer           = void *;
    //!\brief Tag this class as an input access iterator.
    using iterator_category = std::output_iterator_tag;
    //!\}

    /*!\name Constructors, destructor and assignment.
     * \{
     */
    //!\brief Default constructor.
    constexpr out_file_iterator() = default;
    //!\brief Copy constructor.
    constexpr out_file_iterator(out_file_iterator const &) = default;
    //!\brief Copy construction via assignment.
    constexpr out_file_iterator & operator=(out_file_iterator const &) = default;
    //!\brief Move constructor.
    constexpr out_file_iterator (out_file_iterator &&) = default;
    //!\brief Move assignment.
    constexpr out_file_iterator & operator=(out_file_iterator &&) = default;
    //!\brief Use default deconstructor.
    ~out_file_iterator() = default;

    //!\brief Construct with reference to host.
    constexpr out_file_iterator(file_type & _host) noexcept :
        host{&_host}
    {}
    //!\}

    /*!\name Iterator operations
     * \{
     */
    //!\brief This is a no-op, returns reference to self.
    out_file_iterator & operator++()
    {
        return *this;
    }

    //!\brief This is a no-op, returns copy of self. In contrast to input iterators, the return type is required.
    out_file_iterator operator++(int)
    {
        return *this;
    }

    //!\brief Return reference to self.

    out_file_iterator & operator*() noexcept
    {
        return *this;
    }

    /*!\brief Insert the given value into the file, via the file's `push_back()` member.
     * \tparam arg_t    The argument type (the file's `push_back()` will enforce certain constraints).
     * \param[in] arg   The argument to be inserted.
     */
    template <typename arg_t>
    out_file_iterator & operator=(arg_t && arg)
    {
        assert(host != nullptr);
        host->push_back(std::forward<arg_t>(arg));
        return *this;
    }
    //!\}

    /*!\name Comparison operators
     * \brief This iterator never equals its sentinel.
     * \{
     */

    //!\brief Checks whether `*this` is equal to the sentinel (always false).
    constexpr bool operator==(std::ranges::default_sentinel_t const &) const noexcept
    {
        return false;
    }

    //!\brief Checks whether `*this` is not equal to the sentinel (always true).
    constexpr bool operator!=(std::ranges::default_sentinel_t const &) const noexcept
    {
        return true;
    }

    //!\brief Checks whether `it` is equal to the sentinel.
    constexpr friend bool operator==(std::ranges::default_sentinel_t const &,
                                     out_file_iterator const & it) noexcept
    {
        return (it == std::ranges::default_sentinel);
    }

    //!\brief Checks whether `it` is not equal to the sentinel.
    constexpr friend bool operator!=(std::ranges::default_sentinel_t const &,
                                     out_file_iterator const & it) noexcept
    {
        return (it != std::ranges::default_sentinel);
    }
    //!\}

private:
    //!\brief Pointer to file host.
    file_type * host{};
};

} // namespace seqan3::detail
