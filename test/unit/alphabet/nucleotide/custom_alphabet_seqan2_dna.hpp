// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#pragma once

#include <seqan3/alphabet/concept.hpp>

#include <seqan/basic.h>
#include <seqan/modifier/modifier_functors.h>

// Ranks and letters of `seqan2::Dna` and `seqan3::dna4` are the same:
// 'A' = 0, 'C' = 1, 'G' = 2, 'T' = 3.
template <>
struct seqan3::custom::alphabet<seqan2::Dna>
{
    using alphabet_t = seqan2::Dna;

    static constexpr size_t alphabet_size = 4u;

    static uint8_t to_rank(alphabet_t const a) noexcept
    {
        return seqan2::ordValue(a);
    }

    static alphabet_t & assign_rank_to(uint8_t const r, alphabet_t & a) noexcept
    {
        seqan2::assign(a, r);
        return a;
    }

    static char to_char(alphabet_t const a) noexcept
    {
        char c;
        seqan2::assign(c, a);
        return c;
    }

    static alphabet_t & assign_char_to(char const c, alphabet_t & a) noexcept
    {
        seqan2::assign(a, c);
        return a;
    }

    static alphabet_t complement(alphabet_t const a) noexcept
    {
        static seqan2::FunctorComplement<alphabet_t> func;
        return func(a);
    }
};

// Disable `seqan2::begin` overload for `seqan::Dna`.
// Otherwise, `seqan2::Dna` would be considered a range.
// `seqan2::Dna` is indeed not a range, and assuming otherwise breaks many views.
namespace seqan2
{

void begin(seqan2::Dna) = delete;

} // namespace seqan2

// This is only needed if you want to use `seqan3::debug_stream` with `seqan2::Dna`.
namespace seqan3
{

inline debug_stream_type<char> & operator<<(debug_stream_type<char> & stream, seqan2::Dna const & data)
{
    stream << seqan2::convert<char>(data);
    return stream;
}

} // namespace seqan3
