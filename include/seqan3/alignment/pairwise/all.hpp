// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Meta-header for the \link alignment_pairwise Alignment / Pairwise submodule \endlink.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

/*!\defgroup alignment_pairwise Pairwise Alignments
 * \ingroup alignment
 * \brief Provides the algorithmic components for the computation of pairwise alignments.
 *
 * \details
 *
 * # Introduction to pairwise alignment
 *
 * Pairwise sequence alignments can be computed with the free function seqan3::align_pairwise. This function is called
 * in the default case with a sequence pair and an alignment configuration object.
 * Note the type of the pair must model seqan3::tuple_like, e.g. std::tuple or std::pair, with exactly two elements.
 * The algorithm accesses the sequences via the corresponding get interface. Furthermore, the sequences stored in the
 * pair must model std::ranges::viewable_range. This means the type returned by the get interface models either
 * std::ranges::view or std::ranges::range and is an lvalue reference. If you don't know yet what a view or a range is,
 * it is recommended to read through the \ref tutorial_ranges "ranges tutorial".
 * The following code snippet demonstrates a simple use of the pairwise alignment interface.
 *
 * \include doc/tutorial/08_pairwise_alignment/pairwise_alignment_first_global.cpp
 *
 * In this snippet a global alignment over two nucleotide sequences using an edit scoring scheme is computed.
 * Here the helper function std::tie is used to pass the two sequences as a tuple to the alignment algorithm.
 * The special interface of std::tie allows to forward the two sequences as lvalue references such that no copy of the
 * data is involved.
 *
 * There are a lot of applications that need to compute many pairwise sequence alignments. Accordingly, the
 * seqan3::align_pairwise interface offers an overload for ranges over sequence pairs. The following snippet shows
 * a simple use case.
 *
 * \include test/snippet/alignment/pairwise/align_pairwise_range.cpp
 *
 * \attention In addition to the type requirements above the alignment interface requires that the passed sequences
 *            model std::ranges::random_access_range and std::ranges::sized_range in order to work correctly, e.g. a
 *            std::vector.
 * # Configuring pairwise alignments
 *
 * In SeqAn the alignment algorithm can be configured in many different ways. The core of this configuration are the
 * different configuration elements that select specific features of the algorithm. To allow a maximal flexibility
 * the configuration is separated from the alignment interface. This means the algorithm must be configured before it is
 * invoked. The respective alignment configurations are defined in their own
 * namespace called seqan3::align_cfg. This namespace is used to disambiguate configurations for the
 * alignment algorithm with configurations from other algorithms in SeqAn.
 *
 * To compute a pairwise alignment at least two configuration elements must be provided: The alignment **method** and the
 * **scoring scheme**. Thus, a valid alignment configuration must specify what kind of alignment shall be computed, as
 * it strongly depends on the corresponding context. A default wouldn't make much sense here. For similar
 * reasons the scoring scheme has to be always provided by the user.
 *
 * # Global and local alignments
 *
 * There have been numerous adaptions and modifications of the original global alignment problem to solve similar
 * problems such as the local alignment. Here, the goal is to find a maximal homologue region between two
 * sequences that has been conserved during the evolution. Other examples are the **semi-global** alignment which is
 * frequently used in read mapping in order to align a smaller sequence into the context of a larger reference sequence.
 *
 * The standard global and local alignments can be configured using seqan3::align_cfg::method_global and
 * seqan3::align_cfg::method_local, respectively.
 *
 * A **semi-global** alignment can be computed by specifying the free end-gaps in the constructor of
 * seqan3::align_cfg::method_global. The parameters enable, respectively disable, the scoring of
 * leading and trailing gaps at the respective sequence ends (
 * \ref seqan3::align_cfg::free_end_gaps_sequence1_leading "first sequence leading",
 * \ref seqan3::align_cfg::free_end_gaps_sequence2_leading "second sequence leading" or
 * \ref seqan3::align_cfg::free_end_gaps_sequence1_trailing "first sequence trailing",
 * \ref seqan3::align_cfg::free_end_gaps_sequence2_trailing "second sequence trailing" gaps).
 * The SeqAn alignment algorithm allows any free end-gap settings making it a very versatile
 * algorithm.
 * This option, however, is not available for the local alignment where penalising gaps at the ends of the
 * sequences is always disabled.
 *
 * \include{doc} doc/fragments/alignment_configuration_align_config_method_global.md
 *
 * \include{doc} doc/fragments/alignment_configuration_align_config_method_local.md
 *
 * # Using scoring and gap schemes
 *
 * To compute an alignment a scoring and a gap scheme must be provided which give a "score" for substituting, inserting,
 * or deleting a base within the alignment computation. Throughout SeqAn, a positive score implies
 * higher similarity and/or a closer relatedness and a lower or even negative score implies distance.
 * If you are used to dealing with "penalties" or "distances", instead think of "negative scores" when using SeqAn
 * interfaces.
 *
 * ## Scoring two letters
 *
 * Scoring two letters of a single alphabet (or two similar alphabets) is performed by scoring schemes. A scoring
 * scheme is any type that models seqan3::scoring_scheme_for, i.e. it must provide a member function that
 * takes the two letters and returns the scheme-specific score. Algorithms that expect a scoring scheme should check
 * this concept with their respective alphabet(s).
 *
 * Two generic scoring schemes are available in SeqAn:
 *
 *   1. seqan3::nucleotide_scoring_scheme that accepts all nucleotides (and any alphabet that is explicitly
 * convertible to seqan3::dna15)
 *   2. seqan3::aminoacid_scoring_scheme that accepts all amino acids (and any alphabet that is explicitly convertible
 * to seqan3::aa27).
 *
 * These also support scoring two nucleotides/amino acids of different types and they also support modification of
 * their scores via `set_()` functions and by returning references to their internal score matrix. You can however
 * add completely different types, as long as they model seqan3::scoring_scheme_for for the respective sequence
 * alphabet type. In fact, when invoking the seqan3::align_pairwise interface it will be checked at compile time if
 * the provided scoring scheme can be used in combination with the passed sequences and if not a static
 * assertion is raised.
 *
 * The scoring scheme can be configured with the seqan3::align_cfg::scoring_scheme element. Since the scoring scheme is
 * strongly coupled on the sequences to be aligned, there is no default for it. Thus, it is mandatory for
 * the developer to specify this configuration.
 *
 * ## Scoring gaps
 *
 * Throughout SeqAn we use the term gap to refer to an individual gap (see \ref alphabet_gap) and a gap interval to
 * refer to a stretch of consecutive gaps.
 * When aligning two sequences a gap is introduced to mark an insertion or deletion with respect to the other sequence.
 * However, because it is widely recognised that the likelihood of `n` consecutive gaps is much higher than that
 * of `n` individual gaps the scoring of an individual gap or a stretch of gaps is not handled by the scoring scheme.
 * This is based on the assumption that one biological event often introduces more than one gap at a time and
 * that single character gaps are not common due to other biological factors like frame preservation in protein-coding
 * sequences.
 *
 * SeqAn offers the additional seqan3::gap_cost_affine which can be used to set the scores for opening
 * (seqan3::align_cfg::open_score) or extending a gap (seqan3::align_cfg::extension_score).
 *
 * # Combining configuration elements
 *
 * Configurations can be combined using the `|`-operator. If a combination is invalid, a static assertion is raised
 * during the compilation of the program. It will inform the user that some configurations cannot be combined together
 * into one alignment configuration. In general, the same configuration element cannot occur more than once inside of
 * a configuration specification. The following table shows which combinations are possible.
 *
 * | **Config**                                                                  | **0** | **1** | **2** | **3** | **4** | **5** | **6** | **7** | **8** | **9** | **10** | **11** | **12** | **13** | **14** |
 * |:----------------------------------------------------------------------------|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:------:|:------:|:------:|:------:|:------:|
 * | \ref seqan3::align_cfg::band_fixed_size "0: Band"                           |  ❌   |   ✅   |  ✅   |  ✅   |  ✅   |   ✅   |  ✅   |  ✅   |   ✅   |  ✅   |    ✅   |   ✅   |   ✅   |    ✅   |   ✅   |
 * | \ref seqan3::align_cfg::gap_cost_affine "1: Gap scheme affine"              |  ✅   |   ❌   |  ✅   |  ✅   |  ✅   |   ✅   |  ✅   |  ✅   |   ✅   |  ✅   |    ✅   |   ✅   |   ✅   |    ✅   |   ✅   |
 * | \ref seqan3::align_cfg::min_score "2: Min score"                            |  ✅   |   ✅   |  ❌   |  ✅   |  ❌   |   ✅   |  ✅   |  ✅   |   ✅   |  ✅   |    ✅   |   ✅   |   ✅   |    ✅   |   ✅   |
 * | \ref seqan3::align_cfg::method_global "3: Method global"                    |  ✅   |   ✅   |  ✅   |  ❌   |  ❌   |   ✅   |  ✅   |  ✅   |   ✅   |  ✅   |    ✅   |   ✅   |   ✅   |    ✅   |   ✅   |
 * | \ref seqan3::align_cfg::method_local "4: Method local"                      |  ✅   |   ✅   |  ❌   |  ❌   |  ❌   |   ✅   |  ✅   |  ✅   |   ✅   |  ✅   |    ✅   |   ✅   |   ✅   |    ✅   |   ✅   |
 * | \ref seqan3::align_cfg::output_alignment "5: Alignment output"              |  ✅   |   ✅   |  ✅   |  ✅   |  ✅   |   ❌   |  ✅   |  ✅   |   ✅   |  ✅   |    ✅   |   ✅   |   ✅   |    ✅   |   ✅   |
 * | \ref seqan3::align_cfg::output_end_position "6: End positions output"       |  ✅   |   ✅   |  ✅   |  ✅   |  ✅   |   ✅   |  ❌   |  ✅   |   ✅   |  ✅   |    ✅   |   ✅   |   ✅   |    ✅   |   ✅   |
 * | \ref seqan3::align_cfg::output_begin_position "7: Begin positions output"   |  ✅   |   ✅   |  ✅   |  ✅   |  ✅   |   ✅   |  ✅   |  ❌   |   ✅   |  ✅   |    ✅   |   ✅   |   ✅   |    ✅   |   ✅   |
 * | \ref seqan3::align_cfg::output_score "8: Score output"                      |  ✅   |   ✅   |  ✅   |  ✅   |  ✅   |   ✅   |  ✅   |  ✅   |   ❌   |  ✅   |    ✅   |   ✅   |   ✅   |    ✅   |   ✅   |
 * | \ref seqan3::align_cfg::output_sequence1_id "9: Sequence1 id output"        |  ✅   |   ✅   |  ✅   |  ✅   |  ✅   |   ✅   |  ✅   |  ✅   |   ✅   |  ❌   |    ✅   |   ✅   |   ✅   |    ✅   |   ✅   |
 * | \ref seqan3::align_cfg::output_sequence2_id "10: Sequence2 id output"       |  ✅   |   ✅   |  ✅   |  ✅   |  ✅   |   ✅   |  ✅   |  ✅   |   ✅   |  ✅   |    ❌   |   ✅   |   ✅   |    ✅   |   ✅   |
 * | \ref seqan3::align_cfg::parallel "11: Parallel"                             |  ✅   |   ✅   |  ✅   |  ✅   |  ✅   |   ✅   |  ✅   |  ✅   |   ✅   |  ✅   |    ✅   |   ❌   |   ✅   |    ✅   |   ✅   |
 * | \ref seqan3::align_cfg::score_type "12: Score type"                         |  ✅   |   ✅   |  ✅   |  ✅   |  ✅   |   ✅   |  ✅   |  ✅   |   ✅   |  ✅   |    ✅   |   ✅   |   ❌   |    ✅   |   ✅   |
 * | \ref seqan3::align_cfg::scoring_scheme "13: Scoring scheme"                 |  ✅   |   ✅   |  ✅   |  ✅   |  ✅   |   ✅   |  ✅   |  ✅   |   ✅   |  ✅   |    ✅   |   ✅   |   ✅   |    ❌   |   ✅   |
 * | \ref seqan3::align_cfg::vectorised "14: Vectorised"                         |  ✅   |   ✅   |  ✅   |  ✅   |  ✅   |   ✅   |  ✅   |  ✅   |   ✅   |  ✅   |    ✅   |   ✅   |   ✅   |    ✅   |   ❌   |
 *
 * \if DEV
 * There is an additional configuration element \ref seqan3::align_cfg::detail::debug "Debug", which enables the output
 * of the alignment matrices from the DP algorithm using the returned seqan3::alignment_result. It is compatible with
 * all other configuration elements.
 * \endif
 *
 * # Accessing the alignment results
 * \anchor seqan3_align_cfg_output_configurations
 *
 * The seqan3::align_pairwise interface returns a seqan3::algorithm_result_generator_range. This range is a lazy single
 * pass input range over the computed alignments and the range's element types are seqan3::alignment_result objects.
 * Even if only a single alignment is computed a range will be returned since it could be possible that
 * one alignment invocation produces multiple results, e.g. to receive suboptimal alignments.
 * The seqan3::alignment_result object contains only the information that has been requested via the `output`
 * configuration (see below). The algorithm will then choose the most efficient implementation to
 * compute the requested outputs.
 * The following table shows the different output configurations:
 *
 * | **Output option**                                                                        | **Available result**                     |
 * | -----------------------------------------------------------------------------------------|------------------------------------------|
 * | \ref seqan3::align_cfg::output_score "seqan3::align_cfg::output_score"                   | alignment score                          |
 * | \ref seqan3::align_cfg::output_end_position "seqan3::align_cfg::output_end_position"     | end positions of the aligned sequences   |
 * | \ref seqan3::align_cfg::output_begin_position "seqan3::align_cfg::output_begin_position" | begin positions of the aligned sequences |
 * | \ref seqan3::align_cfg::output_alignment "seqan3::align_cfg::output_alignment"           | alignment of the two sequences           |
 * | \ref seqan3::align_cfg::output_sequence1_id "seqan3::align_cfg::output_sequence1_id"     | id of the first sequence                 |
 * | \ref seqan3::align_cfg::output_sequence2_id "seqan3::align_cfg::output_sequence2_id"     | id of the second sequence                |
 *
 * The begin and end positions refer to the begin and end positions of the slices of the original sequences that are
 * aligned. For example, the positions reported for the global alignment correspond to the positions
 * of the original sequences since the entire sequences are encompassed by the global alignment.
 * In case of a local alignment the aligned part might only encompass a part of the original sequences.
 * In this case, the begin and end positions denote the begin and end of the slices of the original sequences that are
 * aligned.
 * To obtain the actual alignment the option seqan3::align_cfg::output_alignment has to be specified.
 * The options can be combiend with each other in order to customise the alignment algorithm and the respective output
 * of the alignment. For example computing the alignment will always incur some run time penalty compared to just
 * computing the score.
 * See the following snippet for some examples:
 *
 * \include test/snippet/alignment/configuration/align_cfg_output_examples.cpp
 *
 * If none of the above configuration was set by the user, then all output options will be enabled by default, i.e.
 * the alignment algorithm will compute every output. Otherwise, if any of the output configurations was set
 * by the user, then only the configured ones are available in the final seqan3::alignment_result.
 * Trying to access an output which has not been configured will raise a static assertion
 * informing the developer about the invalid access.
 *
 * \note Currently, the sequence ids are represented by an internal mechanism and might not refer to the actual id
 *       of the underlying sequences in the respective alignment, rather it is an ongoing number identifying the
 *       computed pair of sequences. In the future, there will be a mechanism for the user to specify the id of
 *       the sequences.
 *
 * # Algorithmic details
 *
 * Since both algorithms ([Smith-Waterman](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm) and
 * [Needleman-Wunsch](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm) algorithm) are based on dynamic
 * programming, they run in quadratic time and memory **O(nm)** (where `n` and `m` are the lengths of the respective
 * aligned sequences).
 *
 * To reduce the time complexity you can use a \ref seqan3::align_cfg::band_fixed_size "banded" alignment. It reduces
 * the runtime by a constant although remaining quadratic and limiting the possible solutions slightly.
 * You can speed up the computation significantly if you \ref seqan3::align_cfg::parallel "parallelize" and simdify your
 * alignment. More about banded and parallelization can be read below.
 *
 * By default a generic alignment algorithm is used that supports all valid alignment configurations but for some
 * special combinations of parameters a notably faster algorithm is available.
 * It is automatically selected if all of the following requirements are satisfied:
 *  * Edit distance gaps, i.e. default initialised seqan3::align_cfg::gap_cost_affine gap scheme.
 *  * Edit distance scoring for \ref alphabet_nucleotide "nucleotide alphabets", i.e. seqan3::align_cfg::scoring_scheme is
 *   initialised with default initialised seqan3::nucleotide_scoring_scheme.
 *  * Global alignment, i.e. seqan3::align_cfg::method_global.
 *
 * There is a special shortcut for the above required scoring/gap configurations called seqan3::align_cfg::edit_scheme,
 * which can be used to safe some typing.
 *
 * The edit configuration can be further specialised with following configs:
 *  * Set a minimal score, i.e. specify the seqan3::align_cfg::min_score configuration
 *  * Compute a semi-global alignment where the end-gaps of the first sequence are free.
 *
 * \note If there was a configuration that is not suitable for the edit distance algorithm the standard alignment
 *       algorithm is executed as a fallback.
 *
 * # Computing banded alignments
 *
 * \include{doc} doc/fragments/alignment_configuration_align_config_band.md
 *
 * # Parallel alignment execution
 *
 * \include{doc} doc/fragments/alignment_configuration_align_config_parallel.md
 *
 * # User callback
 *
 * In some cases, for example when executing the alignments in parallel, it can be beneficial for the performance to
 * use a continuation interface rather than collecting the results first through the seqan3::algorithm_result_generator_range.
 * To be more precise, if more work needs to be done after the alignment has been computed, it could be better to
 * stay within the thread and continue the work rather than buffering the result and computing the next alignment.
 * The alignment algorithm allows the user to specify their own callback function which will be invoked by the alignment
 * algorithm when a seqan3::alignment_result has been computed. To do so, the seqan3::align_cfg::on_result configuration
 * element can be used during the alignment configuration. Note that if seqan3::align_cfg::on_result is specified, the
 * algorithm seqan3::align_pairwise does not return a seqan3::algorithm_result_generator_range anymore. In fact, the algorithm's return
 * type is `void`. The following code snippet illustrates this behavior:
 *
 * \include test/snippet/alignment/pairwise/parallel_align_pairwise_with_callback.cpp
 *
 * \see
 *  - [lecture script - pairwise alignment](https://www.mi.fu-berlin.de/en/inf/groups/abi/teaching/lectures/lectures_past/WS0910/V___Algorithmen_und_Datenstrukturen/scripts/alignment.pdf)\n
 *  - alignment
 */

#pragma once

#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/pairwise/align_result_selector.hpp>
#include <seqan3/alignment/pairwise/alignment_algorithm.hpp>
#include <seqan3/alignment/pairwise/alignment_configurator.hpp>
#include <seqan3/alignment/pairwise/alignment_result.hpp>
#include <seqan3/alignment/pairwise/edit_distance_algorithm.hpp>
#include <seqan3/alignment/pairwise/edit_distance_fwd.hpp>
#include <seqan3/alignment/pairwise/edit_distance_unbanded.hpp>
#include <seqan3/alignment/pairwise/policy/all.hpp>
