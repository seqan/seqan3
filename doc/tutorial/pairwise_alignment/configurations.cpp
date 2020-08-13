//! [include]
#include <seqan3/alignment/configuration/all.hpp>
//! [include]

//! [include_scoring_scheme]
#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
//! [include_scoring_scheme]

//! [include_method]
#include <seqan3/alignment/scoring/method.hpp>
//! [include_method]

//! [include_gap_scheme]
#include <seqan3/alignment/scoring/gap_scheme.hpp>
//! [include_gap_scheme]

//! [include_result]
#include <seqan3/alignment/configuration/align_config_result.hpp>
//! [include_result]

//! [include_band]
#include <seqan3/alignment/configuration/align_config_band.hpp>
//! [include_band]

//! [include_edit]
#include <seqan3/alignment/configuration/align_config_edit.hpp>
//! [include_edit]

#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/alphabet/aminoacid/all.hpp>

int main()
{
{
//! [method_global_free_end_gaps]
    // Example of a semi-global alignment where leading and trailing gaps in the
    // second sequence are not penelised:
    auto config = seqan3::align_cfg::method_global{
                      seqan3::align_cfg::free_end_gaps_sequence1_leading{false},
                      seqan3::align_cfg::free_end_gaps_sequence2_leading{true},
                      seqan3::align_cfg::free_end_gaps_sequence1_trailing{false},
                      seqan3::align_cfg::free_end_gaps_sequence2_trailing{true}};
//! [method_global_free_end_gaps]
(void) config;
}

{
//! [scoring_scheme]
using seqan3::operator""_dna4;
using seqan3::operator""_aa27;

// Define a simple scoring scheme with match and mismatch cost and get the score.
seqan3::nucleotide_scoring_scheme nc_scheme{seqan3::match_score{4}, seqan3::mismatch_score{-5}};
auto sc_nc = nc_scheme.score('A'_dna4, 'C'_dna4); // sc_nc == -5.

// Define a amino acid similarity matrix and get the score.
seqan3::aminoacid_scoring_scheme aa_scheme{};
aa_scheme.set_similarity_matrix(seqan3::aminoacid_similarity_matrix::BLOSUM30);
auto sc_aa = aa_scheme.score('M'_aa27, 'K'_aa27); // sc_aa == 2.
//! [scoring_scheme]
(void) sc_nc;
(void) sc_aa;
}

{
//! [gap_scheme]

// Define a gap scheme with custom gap scores.
seqan3::gap_scheme g{seqan3::gap_score{-1}, seqan3::gap_open_score{-10}};

auto gap = g.get_gap_score();  // gap == -1
auto gap_open = g.get_gap_open_score(); // gap_open == -10
//! [gap_scheme]
(void) gap;
(void) gap_open;
}

{
//! [result]

// Configure the alignment to only compute the score.
auto cfg = seqan3::align_cfg::result{seqan3::with_score};
//! [result]
(void) cfg;
}

{
//! [band]

// Configure a banded alignment.
auto cfg = seqan3::align_cfg::band_fixed_size{seqan3::align_cfg::lower_diagonal{-4},
                                              seqan3::align_cfg::upper_diagonal{4}};
//! [band]
(void) cfg;
}

{
//! [edit]

// Configure an edit distance alignment.
auto cfg = seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme;
//! [edit]
(void) cfg;
}
}
