//! [include]
#include <seqan3/alignment/configuration/all.hpp>
//! [include]

//! [include_aligned_ends]
#include <seqan3/alignment/configuration/align_config_aligned_ends.hpp>
//! [include_aligned_ends]

//! [include_scoring_scheme]
#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
//! [include_scoring_scheme]

//! [include_gap_cost_affine]
#include <seqan3/alignment/configuration/align_config_gap_cost_affine.hpp>
//! [include_gap_cost_affine]

//! [include_output]
#include <seqan3/alignment/configuration/align_config_output.hpp>
//! [include_output]

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
//! [aligned_ends]

seqan3::front_end_first fef{std::true_type{}};
seqan3::back_end_first bef{std::false_type{}};
seqan3::front_end_second fes{true};
seqan3::back_end_second bes{false};

auto cfg_1 = seqan3::align_cfg::aligned_ends{seqan3::end_gaps{fef, bef, fes, bes}};
auto cfg_2 = seqan3::align_cfg::aligned_ends{seqan3::end_gaps{fef, fes}};
//! [aligned_ends]
(void) cfg_1;
(void) cfg_2;
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
//! [gap_cost_affine]

// Define a gap scheme with custom gap scores.
seqan3::align_cfg::gap_cost_affine affine_scheme{seqan3::align_cfg::open_score{-10},
                                                 seqan3::align_cfg::extension_score{-1}};

int open_score = affine_scheme.open_score; // == -10
int extension_score = affine_scheme.extension_score; // == -1
//! [gap_cost_affine]
(void) open_score;
(void) extension_score;
}

{
//! [output]

// Configure the alignment to only compute the score.
auto cfg = seqan3::align_cfg::output_score{};
//! [output]
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
