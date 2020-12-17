// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Meta-header for the \link alignment alignment module \endlink.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

/*!\defgroup alignment Alignment
 * \brief The alignment module contains concepts, algorithms and classes that are related to the computation of
 *        pairwise and multiple sequence alignments.
 *
 * \details
 *
 * # Introduction
 *
 * An essential step in almost every bioinformatics application or pipeline is to determine the evolutionary distances
 * of two or more biological sequences (genomic or protein sequences). To get this information on base level resolution
 * one needs to align these sequences. During this alignment step a score is computed which estimates how similar
 * the sequences in question are. Moreover, an alignment transcript can be computed which describes the insertions,
 * deletions and substitutions of bases necessary to transform one sequence into another.
 *
 * There have been numerous adaptions and modifications of the original global alignment problem to solve similar
 * problems such as the local alignment. Here, the goal is to find a maximal homologue region between two
 * sequences that has been conserved during the evolution. Other examples are the semi-global alignment which is
 * frequently used in read mapping in order to align a smaller sequence into the context of a larger reference sequence.
 *
 * SeqAn offers a generic multi-purpose alignment library comprising all widely known alignment algorithms as well as
 * many special algorithms. These algorithms are all accessible through an easy to use alignment interface which
 * is described below.
 *
 * # Pairwise alignment
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
 * \include doc/tutorial/pairwise_alignment/pairwise_alignment_first_global.cpp
 *
 * In this snippet a global alignment over two nucleotide sequences is computed. Here the helper function std::tie
 * is used to pass the two sequences as a tuple to the alignment algorithm. The special interface of std::tie allows
 * to forward the two sequences as lvalue references such that no copy of the data is involved.
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
 *
 * ## Configuring pairwise alignments
 *
 * In SeqAn the alignment algorithm can be configured in many different ways. The core of this configuration are the
 * different configuration elements that select specific features of the algorithm. To allow a maximal flexibility
 * the configuration is separated from the alignment interface. This means the algorithm must be configured before it is
 * invoked. The respective alignment configurations are defined in their own
 * namespace called seqan3::align_cfg. This namespace is used to disambiguate configurations for the
 * alignment algorithm with configurations from other algorithms in SeqAn.
 * To compute a pairwise alignment at least two configuration elements must be provided: The alignment method and the
 * scoring scheme. Thus, a valid alignment configuration must specify what kind of alignment shall be computed, as
 * it strongly depends on the corresponding context. A default wouldn't make much sense here. For similar
 * reasons the scoring scheme has to be always provided by the user.
 *
 * ### Combining configuration elements
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
 * ## Accessing the alignment results
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
 * ## Using scoring and gap schemes
 *
 * To compute an alignment a scoring and a gap scheme must be provided which give a "score" for substituting, inserting,
 * or deleting a base within the alignment computation. Throughout SeqAn, a positive score implies
 * higher similarity and/or a closer relatedness and a lower or even negative score implies distance.
 * If you are used to dealing with "penalties" or "distances", instead think of "negative scores" when using SeqAn
 * interfaces.
 *
 * ### Scoring two letters
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
 * ### Scoring gaps
 *
 * Throughout SeqAn we use the term gap to refer to an individual gap (see \ref gap) and a gap interval to refer
 * to a stretch of consecutive gaps.
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
 * ## Computing banded alignments
 *
 * SeqAn offers the computation of banded alignments to reduce the running time of the algorithm. This can be
 * helpful if the region in which the optimal alignment exists is known a priori. To specify the banded alignment
 * the developer can use the seqan3::align_cfg::band_fixed_size option.
 * This band configuration is initialised with a seqan3::align_cfg::lower_diagonal and a
 * seqan3::align_cfg::upper_diagonal. The term diagonal is used to describe the position of the band boundary within the
 * alignment matrix. The given value represents the offset that the lower, respectively upper, diagonal is shifted from the
 * main diagonal, which starts in the origin of the alignment matrix. Accordingly, a negative value shifts the band
 * boundary downwards in the alignment matrix and a positive value shifts the band boundary to the right.
 *
 * The band parameters might be restricted depending on the configured alignment algorithm, e.g. the origin of the
 * alignment matrix and the sink (the last cell in the last column) must be covered by the band when a global alignment
 * is ought to be computed.
 * In general, the upper diagonal must always be greater than or equal to the lower diagonal to specify a valid band.
 *
 * ## Global and local alignments
 *
 * The standard global and local alignments can be configured using seqan3::align_cfg::method_global and
 * seqan3::align_cfg::method_local, respectively.
 * A semi-global alignment can be computed by specifying the free end-gaps in the constructor of
 * seqan3::align_cfg::method_global. The parameters enable, respectively disable, the scoring of
 * leading and trailing gaps at the respective sequence ends.
 * The SeqAn alignment algorithm allows any free end-gap settings making it a very versatile
 * algorithm.
 * This option, however, is not available for the local alignment where penalising gaps at the ends of the
 * sequences is always disabled.
 *
 * ## Algorithmic details
 *
 * By default a generic alignment algorithm is used that supports all valid alignment configurations but for some
 * special combinations of parameters a notably faster algorithm is available.
 * It is automatically selected if all of the following requirements are satisfied:
 *  * Edit distance gaps, i.e. default initialised seqan3::align_cfg::gap_cost_affine gap scheme.
 *  * Edit distance scoring for \ref nucleotide "nucleotide alphabets", i.e. seqan3::align_cfg::scoring_scheme is
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
 * # Parallel alignment execution
 *
 * SeqAn's alignment algorithm is internally accelerated using multi-threading. The parallel execution can be selected
 * by specifying the seqan3::align_cfg::parallel configuration element. This will enable the asynchronous execution
 * of the alignments in the backend. For the user interface nothing changes as the returned seqan3::algorithm_result_generator_range
 * will preserve the order of the computed alignment results, i.e. the first result corresponds to the first alignment
 * given by the input range. By default, a thread pool with std::thread::hardware_concurrency many threads will be
 * created on a call to seqan3::align_pairwise and destructed when all alignments have been processed and the
 * seqan3::algorithm_result_generator_range goes out of scope. The configuration element seqan3::align_cfg::parallel can be initialised
 * with a custom thread count which determines the number of threads that will be spawned in the background.
 *
 * ## User callback
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
 */

 #pragma once

 #include <seqan3/alignment/aligned_sequence/all.hpp>
 #include <seqan3/alignment/configuration/all.hpp>
 #include <seqan3/alignment/exception.hpp>
 #include <seqan3/alignment/matrix/all.hpp>
 #include <seqan3/alignment/pairwise/all.hpp>
 #include <seqan3/alignment/scoring/all.hpp>
