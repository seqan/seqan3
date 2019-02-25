// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

//! [create]
#include <seqan3/alphabet/all.hpp> // for working with alphabets directly

using namespace seqan3;

int main ()
{
    // Two characters of seqan3::dna4 alphabet.
    dna4 ade = 'A'_dna4;
    dna4 gua = 'G'_dna4;

    // Two additional characters assigned explicitly from char or rank.
    dna4 cyt, thy;
    assign_char(cyt, 'C');
    assign_rank(thy, 3);

    // Further code here...
//! [create]
    assert(cyt == 'C'_dna4);
    assert(thy == 'T'_dna4);

//! [rank]
    // Get the rank type of an alphabet.
    using rank_type = underlying_rank_t<dna4>;

    // Retrieve the numerical representation (rank) of the characters.
    rank_type rank_a = ade.to_rank();   // => 0
    rank_type rank_g = gua.to_rank();   // => 2
//! [rank]
    assert(rank_a == 0u);
    assert(rank_g == 2u);

//! [size]
    // Get the alphabet size as class member of the alphabet.
    uint8_t const size1 = dna4::value_size;        // => 4

    // Get the alphabet size as global function.
    uint8_t const size2 = alphabet_size_v<dna4>;   // => 4
//! [size]
    assert(size1 == 4u);
    assert(size2 == 4u);

//! [gapped]
    // Assign a gap symbol to a gapped RNA alphabet.
    gapped<rna5> sym = gap{};                         // => -

    // Each seqan3::rna5 character is still valid.
    sym = 'U'_rna5;                                   // => U

    // The alphabet size becomes six (AUGCN-).
    uint8_t const size3 = gapped<rna5>::value_size;   // => 6
//! [gapped]
    assert(size3 == 6u);

    // Prevent -Wunused-variable warnings.
    (void)rank_a;
    (void)rank_g;
    (void)size1;
    (void)size2;
    (void)size3;

//! [closing]
    return 0;
}
//! [closing]

template <typename alphabet_type>
    requires NucleotideAlphabet<alphabet_type>
alphabet_type complement(alphabet_type const alph)
{
    // ...
}
