// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2016, Knut Reinert, FU Berlin
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
// ==========================================================================

#include <meta/meta.hpp>
#include <range/v3/utility/nullptr_v.hpp>

#include <range/v3/iterator_range.hpp>

#include <seqan3/core/meta.hpp>

#pragma once

// namespace seqan3::detail
// {
//
// template <typename t>
// ranges::v3::iterator_range<t, t>
// chunk_type_helper(t *);
//
// template <typename t>
// using chunk_type_ = meta::_t<decltype(detail::chunk_type_helper(ranges::v3::_nullptr_v<t>()))>;
// }  // namespace seqan3::detail

// namespace seqan3
// {
// // Metafunction chunk_type
// template <typename t>
// struct chunk_type : meta::defer<detail::chunk_type_, std::remove_cv_t<t>>
// {};
//
// template <typename t>
// using chunk_type_t = typename chunk_type<t>::type;
// }  // namespace seqan3

namespace seqan3::detail
{

template <typename iterator_t>
class chunk_decorator
{
private:

    chunk_decorator() = default;

    chunk_decorator(chunk_decorator const & /*other*/) = default;
    chunk_decorator(chunk_decorator && /*other*/) = default;

    chunk_decorator & operator=(chunk_decorator const & /*other*/) = default;
    chunk_decorator & operator=(chunk_decorator && /*other*/) = default;

    ~chunk_decorator() = default;

    // Returns itself as derived type.
    iterator_t & derived() noexcept { return static_cast<iterator_t &>(*this); }
    iterator_t const & derived() const noexcept { return static_cast<iterator_t const &>(*this); }

    friend iterator_t;

public:

    template <typename integral_t = size_t>
    inline void
    next_chunk(integral_t const chunk_size = 1)
    {
        derived().next_chunk_impl(chunk_size);
    }

    template <typename integral_t>
    inline void
    advance_chunk(integral_t const offset)
    {
        derived().advance_chunk_impl(offset);
    }

    inline auto
    get_chunk() const noexcept
    {
        return ranges::v3::iterator_range<decltype(derived().chunk_current()), decltype(derived().chunk_end())>
                    {derived().chunk_current(), derived().chunk_end()};
    }
};

template <typename t>
using chunk_type_ = decltype(std::declval<t>().get_chunk());

}  // namespace seqan3::detail

namespace seqan3
{

// Metafunction chunk_type
template <typename t>
struct chunk_type : meta::defer<detail::chunk_type_, std::remove_cv_t<t>>
{};

template <typename t>
using chunk_type_t = typename chunk_type<t>::type;

// Global function wrappers.
template <typename object_t>
inline auto
get_chunk(object_t const & obj)
{
    return obj.get_chunk();
}

template <typename object_t, typename off_t>
inline void
advance_chunk(object_t & obj, off_t const offset)
{
    obj.advance_chunk(offset);
}

template <typename object_t, typename integral_t>
inline void
next_chunk(object_t & obj, integral_t const offset)
{
    obj.next_chunk(offset);
}
}
