// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <seqan3/alphabet/nucleotide/nucleotide_base.hpp>
#include <seqan3/alphabet/nucleotide/rna4.hpp>

// Conditionally add the seqan3 namespace:
// unit/performance tests: in seqan3 namespace
// cookbook entry: no namespace
#ifndef SEQAN3_USE_NAMESPACE
#    define SEQAN3_USE_NAMESPACE 1
#endif

#if SEQAN3_USE_NAMESPACE
namespace seqan3
{
#endif
/* See discussion here: https://github.com/seqan/seqan3/issues/1970
 *
 * If AVX2 is available, this significantly improves the performance
 * (roughly 5 times faster). But if no AVX2 is available. this decreases the
 * performance (roughly 3 times slower).
 */
//![cookbook]
class simd_dna4 : public seqan3::nucleotide_base<simd_dna4, 256>
{
private:
    using base_t = seqan3::nucleotide_base<simd_dna4, 256>;

    friend base_t;         // nucleotide_base
    friend base_t::base_t; // alphabet_base
    friend seqan3::rna4;

public:
    constexpr simd_dna4() noexcept = default;
    constexpr simd_dna4(simd_dna4 const &) noexcept = default;
    constexpr simd_dna4(simd_dna4 &&) noexcept = default;
    constexpr simd_dna4 & operator=(simd_dna4 const &) noexcept = default;
    constexpr simd_dna4 & operator=(simd_dna4 &&) noexcept = default;
    ~simd_dna4() noexcept = default;

    template <std::same_as<seqan3::rna4> t> // template parameter t to accept incomplete type
    constexpr simd_dna4(t const r) noexcept
    {
        assign_rank(r.to_rank());
    }

    using base_t::assign_rank;
    using base_t::base_t;
    using base_t::to_rank;

    static constexpr uint8_t alphabet_size = 4;

    constexpr simd_dna4 & assign_char(char_type const c) noexcept
    {
        char_type const upper_case_char = c & 0b0101'1111;
        rank_type rank = (upper_case_char == 'T') * 3 + (upper_case_char == 'G') * 2 + (upper_case_char == 'C');
        return assign_rank(rank);
    }

    constexpr char_type to_char() const noexcept
    {
        rank_type const rank = to_rank();
        switch (rank)
        {
        case 0u:
            return 'A';
        case 1u:
            return 'C';
        case 2u:
            return 'G';
        default:
            return 'T';
        }
    }

    constexpr simd_dna4 complement() const noexcept
    {
        rank_type rank{to_rank()};
        rank ^= 0b11;
        simd_dna4 ret{};
        return ret.assign_rank(rank);
    }

    static constexpr bool char_is_valid(char_type const c) noexcept
    {
        char_type const upper_case_char = c & 0b0101'1111;
        return (upper_case_char == 'A') || (upper_case_char == 'C') || (upper_case_char == 'G')
            || (upper_case_char == 'T');
    }
};
//![cookbook]

#if SEQAN3_USE_NAMESPACE
}
#endif
