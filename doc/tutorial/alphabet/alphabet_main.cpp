// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <vector>
#include <set>
#include <tuple>

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
    cyt.assign_char('C');
    thy.assign_rank(3);

    // Further code here...
//! [create]
    assert(cyt == 'C'_dna4);
    assert(thy == 'T'_dna4);

//! [rank]
    // Get the rank type of the alphabet (here uint8_t).
    using rank_type = dna4::rank_type;

    // Retrieve the numerical representation (rank) of the characters.
    rank_type rank_a = ade.to_rank();   // => 0
    rank_type rank_g = gua.to_rank();   // => 2
//! [rank]
    assert(rank_a == 0u);
    assert(rank_g == 2u);

//! [char]
    // Get the char type of the alphabet (here char).
    using char_type = dna4::char_type;

    // Retrieve the character representation.
    char_type char_a = ade.to_char();   // => 'A'
    char_type char_g = gua.to_char();   // => 'G'
//! [char]
    assert(char_a == 'A');
    assert(char_g == 'G');

//! [char_strict]
    // Assign an alphabet character with value check.
    cyt.assign_char_strict('C');

    // thy.assign_char_strict('X'); would throw seqan3::invalid_char_assignment
//! [char_strict]
    assert(cyt == 'C'_dna4);

//! [size]
    // Get the alphabet size as class member of the alphabet.
    uint8_t const size1 = dna4::value_size;        // => 4
//! [size]
    assert(size1 == 4u);

//! [gapped]
    // Assign a gap symbol to a gapped RNA alphabet.
    gapped<rna5> sym = gap{};                         // => -

    // Each seqan3::rna5 character is still valid.
    sym = 'U'_rna5;                                   // => U

    // The alphabet size becomes six (AUGCN-).
    uint8_t const size2 = gapped<rna5>::value_size;   // => 6
//! [gapped]
    assert(size2 == 6u);

//! [containers]
    // Examples of different container types with SeqAn's alphabets.
    std::vector<dna5> dna_sequence{"GATTANAG"_dna5};
    std::pair<gapped<dna4>, gapped<dna4>> alignment_column{gap{}, thy};
    std::set<dna4> pyrimidines{'C'_dna4, 'T'_dna4};
//! [containers]

    // Prevent -Wunused-variable warnings.
    (void)rank_a;
    (void)rank_g;
    (void)char_a;
    (void)char_g;
    (void)size1;
    (void)size2;
    (void)dna_sequence;
    (void)alignment_column;
    (void)pyrimidines;

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
