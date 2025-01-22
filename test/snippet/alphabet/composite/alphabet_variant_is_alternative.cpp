// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <gtest/gtest.h>

#include <seqan3/alphabet/composite/alphabet_variant.hpp>
#include <seqan3/alphabet/gap/gap.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>

int main()
{
    using variant_t = seqan3::alphabet_variant<seqan3::dna5, seqan3::gap>;

    static_assert(variant_t::is_alternative<seqan3::dna5>(), "dna5 is an alternative of variant_t");
    static_assert(!variant_t::is_alternative<seqan3::dna4>(), "dna4 is not an alternative of variant_t");
    static_assert(variant_t::is_alternative<seqan3::gap>(), "gap is an alternative of variant_t");
}
