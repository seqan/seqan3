// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
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
// Author: Rene Rahn <rene.rahn AT fu-berlin.de>
// ============================================================================

#pragma once

#ifndef NDEBUG

#include <random>

#include <seqan3/core/concepts/core.hpp>

namespace seqan3::detail
{

struct test_core_concepts_type_a
{};

struct test_core_concepts_type_b : test_core_concepts_type_a
{
    test_core_concepts_type_b(test_core_concepts_type_b const &) = delete;
    test_core_concepts_type_b(test_core_concepts_type_b &&) = default;

    test_core_concepts_type_b & operator=(test_core_concepts_type_b const &) = delete;
    test_core_concepts_type_b & operator=(test_core_concepts_type_b &&) = default;

    template <typename ...args>
    bool operator()(args &&...);
};

struct test_core_concepts_type_c
{
    test_core_concepts_type_c() = default;
    test_core_concepts_type_c(test_core_concepts_type_b const &)
    {}
    explicit test_core_concepts_type_c(test_core_concepts_type_a const &)
    {}

    template <typename ...args>
    void operator()(args &&... );
};

struct test_core_concepts_type_d: test_core_concepts_type_b
{
    test_core_concepts_type_d() = delete;

    test_core_concepts_type_d(test_core_concepts_type_d const &) = delete;
    test_core_concepts_type_d(test_core_concepts_type_d &&) = delete;

    test_core_concepts_type_d & operator=(test_core_concepts_type_d const &) = delete;
    test_core_concepts_type_d & operator=(test_core_concepts_type_d &&) = delete;
    ~test_core_concepts_type_d() = delete;

    template <typename t1, typename t2>
    void operator()(t1 &&, t2 &&);

    template <typename t>
    bool operator()(t &&, t &&);
};

bool operator==(test_core_concepts_type_a const & , test_core_concepts_type_b const &);
bool operator==(test_core_concepts_type_b const & , test_core_concepts_type_a const &);
bool operator==(test_core_concepts_type_b const & , test_core_concepts_type_b const &);
bool operator==(test_core_concepts_type_d const & , test_core_concepts_type_b const &);
bool operator==(test_core_concepts_type_b const & , test_core_concepts_type_d const &);
bool operator==(test_core_concepts_type_d const & , test_core_concepts_type_d const &);
bool operator==(test_core_concepts_type_c const & , test_core_concepts_type_c const &);

bool operator!=(test_core_concepts_type_a const & , test_core_concepts_type_b const &);
bool operator!=(test_core_concepts_type_b const & , test_core_concepts_type_a const &);
bool operator!=(test_core_concepts_type_b const & , test_core_concepts_type_b const &);
bool operator!=(test_core_concepts_type_d const & , test_core_concepts_type_b const &);
bool operator!=(test_core_concepts_type_b const & , test_core_concepts_type_d const &);
bool operator!=(test_core_concepts_type_d const & , test_core_concepts_type_d const &);
bool operator!=(test_core_concepts_type_c const & , test_core_concepts_type_c const &);

bool operator<(test_core_concepts_type_a const & , test_core_concepts_type_a const &);
bool operator<(test_core_concepts_type_a const & , test_core_concepts_type_b const &);
bool operator<(test_core_concepts_type_b const & , test_core_concepts_type_b const &);
bool operator<(test_core_concepts_type_b const & , test_core_concepts_type_a const &);
bool operator>(test_core_concepts_type_a const & , test_core_concepts_type_a const &);
bool operator>(test_core_concepts_type_a const & , test_core_concepts_type_b const &);
bool operator>(test_core_concepts_type_b const & , test_core_concepts_type_b const &);
bool operator>(test_core_concepts_type_b const & , test_core_concepts_type_a const &);
bool operator<=(test_core_concepts_type_a const & , test_core_concepts_type_a const &);
bool operator<=(test_core_concepts_type_a const & , test_core_concepts_type_b const &);
bool operator<=(test_core_concepts_type_b const & , test_core_concepts_type_b const &);
bool operator<=(test_core_concepts_type_b const & , test_core_concepts_type_a const &);
bool operator>=(test_core_concepts_type_a const & , test_core_concepts_type_a const &);
bool operator>=(test_core_concepts_type_a const & , test_core_concepts_type_b const &);
bool operator>=(test_core_concepts_type_b const & , test_core_concepts_type_b const &);
bool operator>=(test_core_concepts_type_b const & , test_core_concepts_type_a const &);

bool operator<(test_core_concepts_type_d const & , test_core_concepts_type_d const &);
bool operator<(test_core_concepts_type_d const & , test_core_concepts_type_b const &);
bool operator<(test_core_concepts_type_b const & , test_core_concepts_type_d const &);
bool operator>(test_core_concepts_type_d const & , test_core_concepts_type_d const &);
bool operator>(test_core_concepts_type_d const & , test_core_concepts_type_b const &);
bool operator>(test_core_concepts_type_b const & , test_core_concepts_type_d const &);
bool operator<=(test_core_concepts_type_d const & , test_core_concepts_type_d const &);
bool operator<=(test_core_concepts_type_d const & , test_core_concepts_type_b const &);
bool operator<=(test_core_concepts_type_b const & , test_core_concepts_type_d const &);
bool operator>=(test_core_concepts_type_d const & , test_core_concepts_type_d const &);
bool operator>=(test_core_concepts_type_d const & , test_core_concepts_type_b const &);
bool operator>=(test_core_concepts_type_b const & , test_core_concepts_type_d const &);

}  // namespace seqan3::detail

// Check same_concept
static_assert(seqan3::same_concept<int, int, int>);
static_assert(!seqan3::same_concept<int, char, int>);

// Check derived_from_conept
static_assert(seqan3::derived_from_conept<seqan3::detail::test_core_concepts_type_b,
                                          seqan3::detail::test_core_concepts_type_a>);
static_assert(!seqan3::derived_from_conept<seqan3::detail::test_core_concepts_type_a,
                                           seqan3::detail::test_core_concepts_type_b>);

// Check implicitly_convertible_to_concept
static_assert(seqan3::implicitly_convertible_to_concept<seqan3::detail::test_core_concepts_type_b,
                                                        seqan3::detail::test_core_concepts_type_c>);
static_assert(!seqan3::implicitly_convertible_to_concept<seqan3::detail::test_core_concepts_type_c,
                                                         seqan3::detail::test_core_concepts_type_b>);
static_assert(!seqan3::implicitly_convertible_to_concept<seqan3::detail::test_core_concepts_type_a,
                                                         seqan3::detail::test_core_concepts_type_c>);

// Check exlicitly_convertible_to_concept
static_assert(seqan3::explicitly_convertible_to_concept<seqan3::detail::test_core_concepts_type_b,
                                                        seqan3::detail::test_core_concepts_type_c>);
static_assert(!seqan3::explicitly_convertible_to_concept<seqan3::detail::test_core_concepts_type_c,
                                                         seqan3::detail::test_core_concepts_type_b>);
static_assert(seqan3::explicitly_convertible_to_concept<seqan3::detail::test_core_concepts_type_a,
                                                        seqan3::detail::test_core_concepts_type_c>);

// Check convertible_to_concept
static_assert(seqan3::convertible_to_concept<seqan3::detail::test_core_concepts_type_b,
                                             seqan3::detail::test_core_concepts_type_c>);
static_assert(!seqan3::explicitly_convertible_to_concept<seqan3::detail::test_core_concepts_type_c,
                                                         seqan3::detail::test_core_concepts_type_b>);
static_assert(!seqan3::convertible_to_concept<seqan3::detail::test_core_concepts_type_a,
                                              seqan3::detail::test_core_concepts_type_c>);

// Check common_reference_concept
static_assert(seqan3::common_reference_concept<int32_t, int16_t, int8_t>);
static_assert(!seqan3::common_reference_concept<int32_t, int16_t, seqan3::detail::test_core_concepts_type_c>);

// Check common_concept
static_assert(seqan3::common_concept<seqan3::detail::test_core_concepts_type_a,
                                     seqan3::detail::test_core_concepts_type_b>);
static_assert(!seqan3::common_reference_concept<seqan3::detail::test_core_concepts_type_a,
                                                seqan3::detail::test_core_concepts_type_c>);

// Check integral_concept
static_assert(seqan3::integral_concept<int>);
static_assert(!seqan3::integral_concept<float>);

// Check signed_integral_concept
static_assert(seqan3::signed_integral_concept<int>);
static_assert(!seqan3::signed_integral_concept<unsigned>);

// Check unsigned_integral_concept
static_assert(!seqan3::unsigned_integral_concept<int>);
static_assert(seqan3::unsigned_integral_concept<unsigned>);

// Check weakly_equality_comparable_concept
static_assert(seqan3::weakly_equality_comparable_concept<seqan3::detail::test_core_concepts_type_a,
                                                         seqan3::detail::test_core_concepts_type_b>);
static_assert(!seqan3::weakly_equality_comparable_concept<seqan3::detail::test_core_concepts_type_a,
                                                          seqan3::detail::test_core_concepts_type_c>);

// Check equality_comparable_concept
static_assert(!seqan3::equality_comparable_concept<seqan3::detail::test_core_concepts_type_a,
                                                   seqan3::detail::test_core_concepts_type_b>);
static_assert(seqan3::equality_comparable_concept<seqan3::detail::test_core_concepts_type_b,
                                                  seqan3::detail::test_core_concepts_type_d>);

// Check weakly_ordered_concept
static_assert(seqan3::weakly_ordered_concept<seqan3::detail::test_core_concepts_type_a,
                                             seqan3::detail::test_core_concepts_type_b>);
static_assert(!seqan3::weakly_ordered_concept<seqan3::detail::test_core_concepts_type_c,
                                              seqan3::detail::test_core_concepts_type_d>);

// Check weakly_ordered_concept
static_assert(!seqan3::totally_ordered_concept<seqan3::detail::test_core_concepts_type_a,
                                               seqan3::detail::test_core_concepts_type_b>);
static_assert(seqan3::totally_ordered_concept<seqan3::detail::test_core_concepts_type_b,
                                              seqan3::detail::test_core_concepts_type_d>);

// Check destructible_concept
static_assert(seqan3::destructible_concept<seqan3::detail::test_core_concepts_type_a>);
static_assert(!seqan3::destructible_concept<seqan3::detail::test_core_concepts_type_d>);

// Check constructible_concept
static_assert(seqan3::constructible_concept<seqan3::detail::test_core_concepts_type_a>);
static_assert(seqan3::constructible_concept<seqan3::detail::test_core_concepts_type_c,
                                            seqan3::detail::test_core_concepts_type_a>);
static_assert(!seqan3::constructible_concept<seqan3::detail::test_core_concepts_type_c,
                                             seqan3::detail::test_core_concepts_type_a,
                                             seqan3::detail::test_core_concepts_type_b>);

// Check default_constructible_concept
static_assert(seqan3::default_constructible_concept<seqan3::detail::test_core_concepts_type_a>);
static_assert(!seqan3::default_constructible_concept<seqan3::detail::test_core_concepts_type_d>);

// Check move_constructible_concept
static_assert(seqan3::move_constructible_concept<seqan3::detail::test_core_concepts_type_b>);
static_assert(!seqan3::move_constructible_concept<seqan3::detail::test_core_concepts_type_d>);

// Check copy_constructible_concept
static_assert(seqan3::copy_constructible_concept<seqan3::detail::test_core_concepts_type_a>);
static_assert(!seqan3::copy_constructible_concept<seqan3::detail::test_core_concepts_type_b>);

// Check movable_concept
static_assert(seqan3::movable_concept<seqan3::detail::test_core_concepts_type_b>);
static_assert(!seqan3::movable_concept<seqan3::detail::test_core_concepts_type_d>);

// Check copyable_concept
static_assert(seqan3::copyable_concept<seqan3::detail::test_core_concepts_type_a>);
static_assert(!seqan3::copyable_concept<seqan3::detail::test_core_concepts_type_b>);

// Check assignable_concept
static_assert(seqan3::assignable_concept<seqan3::detail::test_core_concepts_type_a&,
                                         seqan3::detail::test_core_concepts_type_a const &>);
static_assert(seqan3::assignable_concept<seqan3::detail::test_core_concepts_type_c&,
                                         seqan3::detail::test_core_concepts_type_b const &>);
static_assert(!seqan3::assignable_concept<seqan3::detail::test_core_concepts_type_a&,
                                          seqan3::detail::test_core_concepts_type_c&>);

// Check swappable_concept
static_assert(seqan3::swappable_concept<seqan3::detail::test_core_concepts_type_a&,
                                        seqan3::detail::test_core_concepts_type_a&>);
static_assert(!seqan3::swappable_concept<seqan3::detail::test_core_concepts_type_b,
                                         seqan3::detail::test_core_concepts_type_c>);

// Check semi_regular_concept
static_assert(seqan3::semi_regular_concept<seqan3::detail::test_core_concepts_type_a>);
static_assert(seqan3::semi_regular_concept<seqan3::detail::test_core_concepts_type_c>);
static_assert(!seqan3::semi_regular_concept<seqan3::detail::test_core_concepts_type_b>);
static_assert(!seqan3::semi_regular_concept<seqan3::detail::test_core_concepts_type_d>);

// Check regular_concept
static_assert(!seqan3::regular_concept<seqan3::detail::test_core_concepts_type_a>);
static_assert(!seqan3::regular_concept<seqan3::detail::test_core_concepts_type_b>);
static_assert(seqan3::regular_concept<seqan3::detail::test_core_concepts_type_c>);
static_assert(!seqan3::regular_concept<seqan3::detail::test_core_concepts_type_d>);

// Check invocable_concept
static_assert(!seqan3::invocable_concept<seqan3::detail::test_core_concepts_type_a, int, double,
                                         seqan3::detail::test_core_concepts_type_b>);
static_assert(seqan3::invocable_concept<std::random_device>);
static_assert(seqan3::invocable_concept<seqan3::detail::test_core_concepts_type_c, int, double,
                                        seqan3::detail::test_core_concepts_type_b>);

// Check regular_invocable_concept
static_assert(!seqan3::regular_invocable_concept<seqan3::detail::test_core_concepts_type_a, int, double,
                                                 seqan3::detail::test_core_concepts_type_b>);
//TODO(rrahn): Should not meet the regular_invocable_concept
//static_assert(!seqan3::regular_invocable_concept<std::random_device>);
static_assert(seqan3::regular_invocable_concept<seqan3::detail::test_core_concepts_type_c, int, double,
                                                seqan3::detail::test_core_concepts_type_b>);

// Check predicate_concept
static_assert(!seqan3::predicate_concept<seqan3::detail::test_core_concepts_type_c, int, double,
                                         seqan3::detail::test_core_concepts_type_b>);
static_assert(seqan3::predicate_concept<seqan3::detail::test_core_concepts_type_b, int, double,
                                        seqan3::detail::test_core_concepts_type_b>);

// Check relation_concept
static_assert(!seqan3::relation_concept<seqan3::detail::test_core_concepts_type_d, int, double>);
static_assert(seqan3::relation_concept<seqan3::detail::test_core_concepts_type_d, int, int>);


#endif // NDEBUG
