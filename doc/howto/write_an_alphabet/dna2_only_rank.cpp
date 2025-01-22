// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

//! [struct]
#include <seqan3/alphabet/concept.hpp> // alphabet concept checks

struct dna2
{
    uint8_t rank{};
};
//! [struct]

//! [alphabet_concept]
static_assert(seqan3::alphabet<dna2> == false); // NOT an alphabet
//! [alphabet_concept]

//! [other_concepts]
static_assert(std::copy_constructible<dna2>);       // ok
static_assert(std::totally_ordered<dna2> == false); // NO comparison operators
static_assert(seqan3::semialphabet<dna2> == false); // NOT a semialphabet
static_assert(seqan3::alphabet<dna2> == false);     // NOT an alphabet
//! [other_concepts]
