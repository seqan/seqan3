// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

//! [struct]
#include <seqan3/alphabet/concept.hpp>                   // alphabet concept checks

struct dna2
{
    uint8_t rank{};
};
//! [struct]

//! [alphabet_concept]
static_assert(seqan3::Alphabet<dna2> == false);          // NOT an Alphabet
//! [alphabet_concept]

//! [other_concepts]
static_assert(std::CopyConstructible<dna2>);             // ok
static_assert(std::StrictTotallyOrdered<dna2> == false); // NO comparison operators
static_assert(seqan3::Semialphabet<dna2> == false);      // NOT a Semialphabet
static_assert(seqan3::Alphabet<dna2> == false);          // NOT an Alphabet
//! [other_concepts]
