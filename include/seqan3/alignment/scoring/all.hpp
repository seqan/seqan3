// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Meta-header for the \link scoring Scoring sub-module \endlink.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

/*!\defgroup scoring Scoring
 * \brief Provides the data structures used for scoring alphabets and sequences.
 * \ingroup alignment
 *
 * \details
 *
 * ### Nomenclature
 *
 * Throughout SeqAn3 we refer to the "score" of two letters or two sequences and imply that a higher score indicates
 * higher similarity and/or a closer relatedness. A lower or even negative score implies distance.
 * If you are used to dealing with "penalties" or "distances", instead think of "negative scores" when using SeqAn3
 * interfaces.
 *
 * ### Scoring two letters
 *
 * Scoring two letters of a single alphabet (or two similar alphabets) is performed by scoring schemes. A scoring
 * scheme is any type that models seqan3::ScoringScheme, i.e. it must provide a member function that
 * takes the two letters and returns the scheme-specific score. Algorithms that expect a scoring scheme should check
 * this concept with their respective alphabet(s).
 *
 * Two generic scoring schemes are provided:
 *
 *   1. seqan3::nucleotide_scoring_scheme that accepts all nucleotides (and any alphabet that is explicitly
 * convertible to seqan3::dna15)
 *   2. seqan3::aminoacid_scoring_scheme that accepts all amino acids (and any alphabet that is explicitly convertible
 * to seqan3::aa27).
 *
 * These also support scoring two nucleotides/amino acids of different types and they also support modification of
 * their scores via `set_()` functions and by returning references to their internal score matrix. You can however
 * add completely different types, as long as they model seqan3::ScoringScheme.
 *
 * <small><i>The base type seqan3::scoring_scheme_base is only implementation detail and not required for most users.
 * Types that model seqan3::ScoringScheme can (but don't need to!) inherit from it.</i></small>
 *
 * ### Scoring gaps
 *
 * Gaps are usually not scored on a per-character-basis, because it is widely recognised that the likelihood of
 * `n` consecutive gaps is much higher than that of `n` individual gaps. This is based on the assumption that
 * one biological event often introduces more than one gap at a time and that single character gaps are not common
 * due to other bioligical factors like frame preservation in protein-coding sequences.
 *
 * Therefore the default scoring schemes do not cover the seqan3::gap alphabet and there is an additional
 * seqan3::gap_scheme whose score() function takes the length as parameter. This forces you to handle gaps separately
 * in your algorithms.
 */

#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/gap_scheme.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/scoring_scheme_base.hpp>
#include <seqan3/alignment/scoring/scoring_scheme_concept.hpp>
