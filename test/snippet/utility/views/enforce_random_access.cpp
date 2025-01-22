// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <string>

#include <seqan3/alignment/decorator/gap_decorator.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/views/enforce_random_access.hpp>

int main()
{
    using namespace seqan3::literals;

    // A gap decorator is a pseudo random access range using logarithmic time complexity internally.
    auto seq = "ACGTACGACT"_dna4;
    seqan3::gap_decorator aligned_seq{seq};

    // It does not fulfill random access semantics because it does not allow constant time access to arbitrary
    // elements in the range. Thus, it is only a bidirectional range by default.
    static_assert(std::ranges::bidirectional_range<decltype(aligned_seq)>);
    static_assert(!std::ranges::random_access_range<decltype(aligned_seq)>);

    // The default interface models at most std::bidirectional_range.
    auto it = std::ranges::begin(aligned_seq); // Returned iterator models std::bidirectional_iterator.
    std::ranges::advance(it, 3);               // Advancing the iterator takes linear time.
    seqan3::debug_stream << *it << '\n';       // "T"

    // After adapting it with enforce_random_access, the returned range models std::random_access_range.
    auto aligned_seq_ra = aligned_seq | seqan3::views::enforce_random_access;

    // The pesudo_random_access wrapper returns a view that enforces random_access.
    // Note that this does not mean that the semantic requirements have been changed. Only all syntactical
    // interfaces for the std::ranges::random_access_range are modeled now by the returned view.
    // Access time still depends on the underlying range.
    static_assert(std::ranges::random_access_range<decltype(aligned_seq_ra)>);

    auto it_ra = std::ranges::begin(aligned_seq_ra); // Returned iterator models std::random_access_iterator.
    std::ranges::advance(it_ra, 3);         // Advancing the iterator takes now logarithmic time instead linear.
    seqan3::debug_stream << *it_ra << '\n'; // "T"
}
