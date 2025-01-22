# Pairwise Alignment {#tutorial_pairwise_alignment}

<!-- SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
     SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
     SPDX-License-Identifier: CC-BY-4.0
-->

<b>Learning Objective:</b> <br/>

In this tutorial, you will learn how to compute pairwise sequence alignments with SeqAn.
This tutorial is a walkthrough with links to the API documentation and is also meant as a source for copy-and-paste code.

\tutorial_head{Intermediate, 60-90 min, \ref setup \ref tutorial_alphabets \ref tutorial_ranges, }

[TOC]

---

# Introduction

Aligning biological sequences is a very prominent component in many bioinformatics applications and pipelines.
Well known genomic applications that use pairwise sequence alignments are read mapping, genome assembly, variant
detection, multiple sequence alignment as well as protein search.

The goal of the pairwise alignment is to obtain an optimal transcript that describes how two DNA sequences are related
to each other by means of substitutions, insertions, or deletions. The computed transcript describes then the operations
necessary to translate the one sequence into the other, as can be seen in the following picture.

\image html align_transcript.png width=800px

The alignment problem is solved with a dynamic programming (DP) algorithm which runs in \f$ (\mathcal{O}(n^2))\f$ time
and space. Besides the global alignment approach, many more variations of this DP based algorithm have been developed
over time. SeqAn unified all of these approaches into a single DP core implementation which can be extended easily and
thus, with all possible configurations, is a very versatile and powerful tool to compute many desired alignment variants.

# Computing pairwise alignments

Let us first have a look at an example of computing a global alignment in SeqAn.

\includelineno doc/tutorial/08_pairwise_alignment/pairwise_alignment_first_global.cpp

In the above example, we want to compute the similarity of two seqan3::dna4 sequences using a global alignment.
For this, we need to import the dna4 alphabet, the seqan3::hamming_scoring_scheme and the seqan3::align_pairwise in
the beginning of our file. We also import the seqan3::debug_stream which allows us to print various types in a nice
formatted manner.

At the beginning of the file, we are defining our two DNA sequences `s1` and `s2`.
If you feel puzzled about what the `_dna4` suffix does, we recommend reading the \ref tutorial_alphabets and
\ref tutorial_ranges before.
In line 16-17 we configure the alignment job with the most simplistic configuration possible.
In this case, it is a global alignment with edit distance, i.e. mismatches, insertions and deletions are scored with `-1`
and matches are scored with `0`.
Later in this tutorial we will give a more detailed description of the \ref alignment_configurations "configuration" and
how it can be used.
The minimum requirement for computing a pairwise alignment is to specify the alignment method
(seqan3::align_cfg::method_local or seqan3::align_cfg::method_global) and the
seqan3::align_cfg::scoring_scheme configuration elements. The first one selects the internal algorithm and the second one
provides the scoring scheme that should be used to score a pair of sequence characters.

Now we are going to call seqan3::align_pairwise. This interface requires two arguments: a tuple or a range of tuples
with exactly two elements and the configuration object. Independent of the number of pairwise alignment jobs submitted
to the algorithm it always returns a range over seqan3::alignment_result. Later we will see how we can use a
continuation interface that calls a user-defined function rather than iterating over the results sequentially.
Finally, we output the score for the computed alignment.

\attention Just calling `pairwise_align` does nothing as it returns a range that is evaluated in a lazy manner.
           Only when calling begin or incrementing the iterator over the range the alignment computation is invoked.

\assignment{Assignment 1}
Copy and paste the minimal working example from above into a cpp file in your working directory and compile and run it.
Add the two sequences "ACGTGACTGACT" and "AGGTACGAGCGACACT" to the set and compute all pairwise sequence alignments
of the four sequences and output the scores.

\hint
You can use the seqan3::views::pairwise_combine to generate all pairwise combinations.
\endhint
\endassignment
\solution

\include doc/tutorial/08_pairwise_alignment/pairwise_alignment_solution_1.cpp

First, we create the vector of seqan3::dna4 sequences. We keep the configuration as is and then modify the initial code
to a range-based for loop looping over the alignment results. Since the seqan3::alignment_result is a class template and the
template parameters are determined during the configuration step we use auto as the result type.
The current result is cached inside of the lazy range and we capture the result as `const &` in order to not tamper with
the result values.

\endsolution

Congratulations, you have computed your first pairwise alignment with SeqAn!
As you can see, the interface is really simple, yet the configuration object makes it extremely flexible to conduct
various different alignment calculations. In the following chapter, you will learn more about the various configuration
possibilities.

# Alignment configurations
\anchor alignment_configurations

The configuration object is the core of the alignment interface. It allows to easily configure the alignment algorithm
without changing the interface of the pairwise alignment. It uses a seqan3::configuration object to chain different
configuration elements together using the logical or-operator ('|'-operator).
You can find an overview over the available configurations \ref seqan3::configuration "here".
The configurations for the alignment module are available in:

\snippet doc/tutorial/08_pairwise_alignment/configurations.cpp include

The configuration elements are all classes that wrap the actual information necessary for the configuration of the
alignment algorithm. Depending on the configuration specification certain features of the algorithm are enabled or
disabled. Moreover, during the initialisation of the algorithm, the best implementation is chosen based on the given
configurations. To avoid possible ambiguities with the configurations of other algorithms, the configuration elements
for the alignment are defined in the special namespace seqan3::align_cfg.

## Global and semi-global alignment

The most straightforward algorithm is the global alignment which can be configured using the
seqan3::align_cfg::method_global.

\note The method configuration must be given by the user as it strongly depends on the application context.
      It would be wrong for us to assume what the intended default behaviour should be.

The global alignment can be further refined by initialising the seqan3::align_cfg::method_global configuration element
with the free end-gap specifiers. They specify whether gaps at the end of the sequences are penalised.
In SeqAn you can configure this behaviour for every end, namely for leading and trailing gaps of the first and second
sequence. seqan3::align_cfg::method_global is constructed with 4 free end-gap specifiers (one for every end):

 - seqan3::align_cfg::free_end_gaps_sequence1_leading - If set to true, aligning leading gaps in first sequence is
                                                        not penalised.
 - seqan3::align_cfg::free_end_gaps_sequence2_leading - If set to true, aligning leading gaps in second sequence is
                                                        not penalised.
 - seqan3::align_cfg::free_end_gaps_sequence1_trailing - If set to true, aligning trailing gaps in first sequence is
                                                        not penalised.
 - seqan3::align_cfg::free_end_gaps_sequence2_trailing - If set to true, aligning trailing gaps in second sequence is
                                                        not penalised.

The following code snippet demonstrates the different use cases:

\snippet doc/tutorial/08_pairwise_alignment/configurations.cpp include_method
\snippet doc/tutorial/08_pairwise_alignment/configurations.cpp method_global_free_end_gaps

The order of arguments is fixed and must always be as shown in the example.

\assignment{Assignment 2}

Adapt the code from Assignment 1 to compute the semi-global alignment for the case where both ends of the first
sequence can be aligned to gaps without penalising them. Note that in such a semi-global alignment, the first sequence
would be aligned as an infix of the second sequence.

\endassignment
\solution

\include doc/tutorial/08_pairwise_alignment/pairwise_alignment_solution_2.cpp

To accomplish our goal we initialise the `method_global` option with the free end specifiers
for sequence 1 set to `true`, and those for sequence 2 with `false`.

\endsolution

## Scoring schemes and gap schemes

A scoring scheme can be queried to get the score for substituting two alphabet values using the `score` member function.
Currently, SeqAn supports three scoring schemes. Besides seqan3::hamming_scoring_scheme, there is the
seqan3::nucleotide_scoring_scheme used for aligning nucleotide sequences and
seqan3::aminoacid_scoring_scheme used for aligning amino acid sequences.
You can import them with the following includes:

\snippet doc/tutorial/08_pairwise_alignment/configurations.cpp include_scoring_scheme

In many biological applications, using the seqan3::hamming_scoring_scheme
([Hamming Distance](https://en.wikipedia.org/wiki/Hamming_distance)) is not enough as it entirely ignores the biological
background of sequence alignments.
If required, you can use the nucleotide scoring scheme when aligning \ref alphabet_nucleotide "nucleotides" and the
the amino acid scoring scheme when working with \ref alphabet_aminoacid "aminoacids".
These scoring schemes allow a finer control of scoring two aligned letters.
For example, you can set different match and mismatch scores for different pairs of aligned letters.
By default, the scoring schemes are initialised by
setting a \ref seqan3::scoring_scheme_base::set_simple_scheme "simple scheme" consisting of seqan3::match_score and
seqan3::mismatch_score.
But it is also possible to provide a \ref seqan3::scoring_scheme_base::set_custom_matrix "custom matrix".
The amino acid scoring scheme can additionally be \ref seqan3::aminoacid_scoring_scheme::set_similarity_matrix
"initialised" with a predefined substitution matrix that can be accessed via the seqan3::aminoacid_similarity_matrix
enumeration class.

\snippet doc/tutorial/08_pairwise_alignment/configurations.cpp scoring_scheme

\note You can also provide your own scoring scheme implementation. It only has to model the seqan3::scoring_scheme_for concept.

Similarly to the scoring scheme, you can use the seqan3::align_cfg::gap_cost_affine to customise the gap penalties used for
the alignment computation. The default initialised seqan3::align_cfg::gap_cost_affine sets the score for a gap to `-1`
and for a gap opening to `0`. Note that the gap open score is added to the gap score when a gap is opened within the
alignment computation. Therefore, setting the gap open score to `0` disables affine gaps altogether.
You can pass a seqan3::align_cfg::extension_score and a seqan3::align_cfg::open_score object to initialise the scheme
with custom gap penalties. The penalties can be changed later by using the respective member variables
`extension_score` and `open_score`.

\attention SeqAn's alignment algorithm computes the maximal similarity score, thus the match score must be set to a
positive value and the score for mismatch and gaps must be negative in order to maximise over the matching letters.

\snippet doc/tutorial/08_pairwise_alignment/configurations.cpp include_gap_cost_affine
\snippet doc/tutorial/08_pairwise_alignment/configurations.cpp gap_cost_affine

To configure the scoring scheme and the gap scheme for the alignment algorithm you need to use the
seqan3::align_cfg::scoring_scheme and the seqan3::align_cfg::gap_cost_affine configurations. The
seqan3::align_cfg::scoring_scheme is mandatory - similarly to the alignment method configuration. It would be
wrong to assume what the default scoring scheme should be. If you do not provide these configurations, the compilation
will fail with a corresponding error message. Not providing the gap scheme is ok. In this case, the default initialised
gap scheme will be used for the alignment computation.

\assignment{Assignment 3}
Compute the alignment of the two amino acid sequences listed below using an affine gap scheme with gap costs of `-2` and
gap open costs of `-9`. Use the BLOSUM62 similarity matrix.
 -  QFSEEILSDIYCWMLQCGQERAV
 -  AFLPGWQEENKLSKIWMKDCGCLW
\endassignment
\solution

\include doc/tutorial/08_pairwise_alignment/pairwise_alignment_solution_3.cpp

Only a few parts of our algorithm need to be adapted. First, we use an amino acid scoring scheme and initialise it with
the respective similarity matrix. Second, we initialise the gap scheme to represent the affine gap model as given in
the assignment. Et voilà, we have computed a pairwise alignment over aminoacid sequences.
\endsolution

## Alignment result

So far, we have only used the score. However, in many situations the final alignment is required, e.g. when
mapping reads and the user wishes to write the alignment to the final SAM/BAM file.
In SeqAn you can simply configure what is going to be computed by the alignment algorithm using the
different \ref seqan3_align_cfg_output_configurations "output configurations".

\snippet doc/tutorial/08_pairwise_alignment/configurations.cpp include_output
\snippet doc/tutorial/08_pairwise_alignment/configurations.cpp output

Accordingly, the alignment algorithm is configured to use the best implementation to obtain the desired result.
The following table shows the different outcomes that can be configured:

| **Output option**                                                                        | **Available result**                     |
| -----------------------------------------------------------------------------------------|------------------------------------------|
| \ref seqan3::align_cfg::output_score "seqan3::align_cfg::output_score"                   | alignment score                          |
| \ref seqan3::align_cfg::output_end_position "seqan3::align_cfg::output_end_position"     | end positions of the aligned sequences   |
| \ref seqan3::align_cfg::output_begin_position "seqan3::align_cfg::output_begin_position" | begin positions of the aligned sequences |
| \ref seqan3::align_cfg::output_alignment "seqan3::align_cfg::output_alignment"           | alignment of the two sequences           |
| \ref seqan3::align_cfg::output_sequence1_id "seqan3::align_cfg::output_sequence1_id"     | id of the first sequence                 |
| \ref seqan3::align_cfg::output_sequence2_id "seqan3::align_cfg::output_sequence2_id"     | id of the second sequence                |

The final result is returned as a seqan3::alignment_result object. This object offers special member functions to access
the stored values. If you try to access a value, e.g. the alignment, although you didn't specify
seqan3::align_cfg::output_alignment in the output configuration, a static assertion will be triggered during
compilation.

\note If you don't specify any of the above mentioned output configurations then by default all options are enabled and
      will be computed. In order to potentially increase the performance of the alignment algorithm only enable those
      options that are needed for your use case.

\assignment{Assignment 4}
Compute the overlap alignment of the following two sequences. Use a linear gap scheme with a gap score of `-4` and
a simple scoring scheme with mismatch `-2` and match `4`.

 -  TTACGTACGGACTAGCTACAACATTACGGACTAC
 -  GGACGACATGACGTACGACTTTACGTACGACTAGC

\hint
An alignment is called overlap alignment if it allows free end-gaps at each sequence end. With this configuration
overlaps between two sequences can be computed which is a common use case during the overlap layout consensus assembly.
\endhint
\endassignment
\solution

\include doc/tutorial/08_pairwise_alignment/pairwise_alignment_solution_4.cpp

\endsolution

## Banded alignment

In many situations it is not necessary to compute the entire alignment matrix but only a part of it. This has
positive impacts on the performance. To limit the computation space the alignment matrix can be bounded by a band.
Thus, only the alignment is computed that fits in this band. Notably, this alignment does not need to be the optimal alignment.
However, in many cases, we can give a rough bound on how similar the sequences will be and therefor use the banded alignment.
To do so, you can configure the alignment using the seqan3::align_cfg::band_fixed_size option. This configuration
element will be initialised with a seqan3::align_cfg::lower_diagonal and seqan3::align_cfg::upper_diagonal parameter.

\snippet doc/tutorial/08_pairwise_alignment/configurations.cpp include_band
\snippet doc/tutorial/08_pairwise_alignment/configurations.cpp band

\assignment{Assignment 5}
Use the example from assignment 4 and compute it in a band with lower diagonal set to `-3` and upper diagonal set to `8`.
How does the result change?

\endassignment
\solution

\include doc/tutorial/08_pairwise_alignment/pairwise_alignment_solution_5.cpp

\endsolution

## Edit distance

A special form of the pairwise sequence alignment is the edit distance. This distance metric counts the number of
edits necessary to transform one sequence into the other. The cost model for the edit distance is fixed. In particular,
the match score is `0` and the scores for a mismatch and a gap is `-1`. Due to the special metric, a fast
[bitvector](https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.332.9395&rep=rep1&type=pdf) implementation can be
used to compute the edit distance. This happens in SeqAn automatically if the respective configurations are used.
To do so, you need to explicitly use the seqan3::hamming_scoring_scheme as well as the default gap scheme
initialised with `-1` for a gap and `0` for a gap open score while computing a seqan3::global_alignment.
To make the configuration easier, we added a shortcut called seqan3::align_cfg::edit_scheme.

\snippet doc/tutorial/08_pairwise_alignment/configurations.cpp include_edit
\snippet doc/tutorial/08_pairwise_alignment/configurations.cpp edit

The `edit_scheme` still has to be combined with an alignment method. When combining it
with the seqan3::align_cfg::method_global configuration element, the edit distance algorithm
can be further refined with free end-gaps (see section `Global and semi-global alignment`).

\attention Only the following free end-gap configurations are supported for the
global alignment configuration with the edit scheme:
- no free end-gaps (all free end-gap specifiers are set to `false`)
- free end-gaps for the first sequence (free end-gaps are set to `true` for the first and
  to `false` for the second sequence)
Using any other free end-gap configuration will disable the edit distance algorithm, i.e. the fast bitvector
algorithm, and will fall back to the standard pairwise alignment.

### Refine edit distance

The edit distance can be further refined using the seqan3::align_cfg::min_score configuration to fix an edit score
(a limit of the allowed number of edits). If the respective alignment could not find a solution within the given error
bound, the resulting score is infinity (corresponds to std::numeric_limits::max). Also the alignment and the begin and
end positions of the alignment can be computed using a combination of the seqan3::align_cfg::output_alignment,
seqan3::align_cfg::output_begin_position and seqan3::align_cfg::output_end_position options.

\assignment{Assignment 6}

Compute all pairwise alignments from the assignment 1 (only the scores). Only allow at most 7 errors and
filter all alignments with 6 or fewer errors.

\hint
You can use the std::views::filter to get only those alignments that fit the requirements.
\endhint

\endassignment
\solution

\include doc/tutorial/08_pairwise_alignment/pairwise_alignment_solution_6.cpp

\endsolution

## Invalid configurations

Chaining the configurations to build an individual alignment algorithm is a strong advantage of this design. However,
some combinations would result in an invalid alignment configuration. To explicitly prevent this we added some security
details. First, if a combination is invalid (for example by providing the same configuration more than once) a static
assert will inform you about the invalid combination. \ref seqan3::configuration "Here" you can find a
table depicting the valid configurations. Further, if the seqan3::align_pairwise is called, it checks if the input
data can be used with the given configuration. For example, a static assertion is emitted if the alphabet types of the
sequences together with the provided scoring scheme do not model the concept seqan3::scoring_scheme_for.
Other possible errors are invalid band settings where the initialised band does not intersect with the actual alignment
matrix (the lower diagonal starts beyond the end of the first sequence).

<!-- # Parallel execution -->
