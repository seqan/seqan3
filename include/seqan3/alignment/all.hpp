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
 * std::ranges::view or std::ranges::range and is an lvalue reference. If you don't know yet what a view or a range is
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
 * In addition to the type requirements above the alignment interface requires std::ranges::random_access_range and
 * std::ranges::sized_range in order to work correctly.
 *
 * ## Configuring pairwise alignments
 *
 * In SeqAn the alignment algorithm can be configured in many different ways. The core of this configuration are the
 * different configuration elements that select specific features of the algorithm. To allow a maximal flexibility
 * the configuration is separated from the alignment interface. This means that before the alignment algorithm
 * is invoked, the algorithm must be configured. The respective alignment configurations are defined in their own
 * namespace called seqan3::align_cfg. This namespace is used to disambiguate configurations for the
 * alignment algorithm with configurations from other algorithms in SeqAn.
 * To compute a pairwise alignment at least two configuration elements must be provided, namely the
 * the alignment method and the seqan3::align_cfg::scoring.
 *
 * ### Combining configuration elements
 *
 * Configurations can be combined using the `|`-operator. If a combination is invalid, a static assertion is triggered
 * during compilation and will inform the user that the the last config cannot be combined with any of the configs from
 * the left-hand side of the configuration specification. Unfortunately, the names of the invalid
 * types cannot be printed within the static assert, but the following table shows which combinations are possible.
 * In general, the same configuration element cannot occur more than once inside of a configuration specification.
 *
 * | **Config**                                               | **0** | **1** | **2** | **3** | **4** | **5** | **6** | **7** | **8** | **9** |
 * |:---------------------------------------------------------|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|
 * | \ref seqan3::align_cfg::aligned_ends "0: Aligned ends"   |  ❌   |  ✅   |  ✅   |  ✅   |  ❌   |  ✅   |  ✅   |  ✅   |  ✅   |  ✅   |
 * | \ref seqan3::align_cfg::band_fixed_size "1: Band"        |  ✅   |  ❌   |  ✅   |  ✅   |  ✅   |  ✅   |  ✅   |  ✅   |  ✅   |  ✅   |
 * | \ref seqan3::align_cfg::gap "2: Gap scheme"              |  ✅   |  ✅   |  ❌   |  ✅   |  ✅   |  ✅   |  ✅   |  ✅   |  ✅   |  ✅   |
 * | \ref seqan3::align_cfg::max_error "3: Max error"         |  ✅   |  ✅   |  ✅   |  ❌   |  ✅   |  ❌   |  ✅   |  ✅   |  ✅   |  ✅   |
 * | \ref seqan3::align_cfg::method_global "4: Method global" |  ❌   |  ✅   |  ✅   |  ✅   |  ❌   |  ❌   |  ✅   |  ✅   |  ✅   |  ✅   |
 * | \ref seqan3::align_cfg::method_local "5: Method local"   |  ✅   |  ✅   |  ✅   |  ❌   |  ❌   |  ❌   |  ✅   |  ✅   |  ✅   |  ✅   |
 * | \ref seqan3::align_cfg::parallel "6: Parallel"           |  ✅   |  ✅   |  ✅   |  ✅   |  ✅   |  ✅   |  ❌   |  ✅   |  ✅   |  ✅   |
 * | \ref seqan3::align_cfg::result "7: Result"               |  ✅   |  ✅   |  ✅   |  ✅   |  ✅   |  ✅   |  ✅   |  ❌   |  ✅   |  ✅   |
 * | \ref seqan3::align_cfg::scoring "8: Scoring scheme"      |  ✅   |  ✅   |  ✅   |  ✅   |  ✅   |  ✅   |  ✅   |  ✅   |  ❌   |  ✅   |
 * | \ref seqan3::align_cfg::vectorised "9: Vectorised"       |  ✅   |  ✅   |  ✅   |  ✅   |  ✅   |  ✅   |  ✅   |  ✅   |  ✅   |  ❌   |
 *
 * \if DEV
 * There is an additional configuration element \ref seqan3::align_cfg::debug "Debug", which enables the output of the
 * alignment matrices from the DP algorithm using the returned seqan3::alignment_result. It is compatible with all
 * other configuration elements.
 * \endif
 *
 * ## Accessing the computed alignment
 *
 * The seqan3::align_pairwise interface returns a seqan3::algorithm_result_generator_range. This range is a single pass range over
 * the computed alignments and the range's element types are seqan3::alignment_result objects.
 * Even if only a single alignment is computed a range will be returned since it could be possible that
 * one alignment invocation produces multiple results, e.g. to receive suboptimal alignments.
 * The seqan3::alignment_result object contains only the information that has been requested via the alignment
 * configuration. The seqan3::align_cfg::result configuration element allows to limit the parts that are computed
 * in the alignment algorithm. The following table shows the parts that are computed depending on the
 * seqan3::align_cfg::result configuration:
 *
 * | **Entity**                                                                   | **Available result**                                  |
 * | -----------------------------------------------------------------------------|------------------------------------------------------ |
 * | \ref seqan3::align_cfg::result::with_score  "with_score"                     | alignment score                                       |
 * | \ref seqan3::align_cfg::result::with_end_positions  "with_end_positions"     | alignment score; back coordinate                      |
 * | \ref seqan3::align_cfg::result::with_begin_positions  "with_begin_positions" | alignment score; back and front coordinate            |
 * | \ref seqan3::align_cfg::result::with_alignment  "with_alignment"             | alignment score; back and front coordinate; alignment |
 *
 * The back coordinate stores the end of the alignment within both sequences. Note that theses positions are
 * inclusive. The front coordinate stores the begin of the alignment part in both sequences.
 * With the last option the complete alignment is computed. If seqan3::align_cfg::result is not specified only the
 * score is available. Accessing one of the results that have not been requested will trigger a static assertion
 * informing the developer about the invalid access.
 *
 * ## Using scoring and gap schemes
 *
 * To compute an alignment a scoring and a gap scheme must be provided which give a "score" for substituting, inserting,
 * or deleting a base within the alignment computation. Throughout SeqAn a positive score implies
 * higher similarity and/or a closer relatedness and a lower or even negative score implies distance.
 * If you are used to dealing with "penalties" or "distances", instead think of "negative scores" when using SeqAn
 * interfaces.
 *
 * ### Scoring two letters
 *
 * Scoring two letters of a single alphabet (or two similar alphabets) is performed by scoring schemes. A scoring
 * scheme is any type that models seqan3::scoring_scheme, i.e. it must provide a member function that
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
 * add completely different types, as long as they model seqan3::scoring_scheme.
 *
 * The scoring scheme can be configured with the seqan3::align_cfg::scoring element. Since the scoring scheme is
 * strongly coupled on the sequences to be aligned it can not be defaulted. Thus, it is mandatory for the
 * the developer to specify the seqan3::align_cfg::scoring configuration.
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
 * SeqAn offers the additional seqan3::gap_scheme which can be used to set the scores for opening or extending a gap.
 *
 * The gap scheme can be configured with the seqan3::align_cfg::gap element. If the configuration is not specified,
 * the algorithm uses edit distance scores (`-1`) for deletion/insertion.
 *
 * ## Computing banded alignments
 *
 * SeqAn offers the computation of banded alignments to reduce the running time of the algorithm. This can be
 * helpful if the region in which the optimal alignment exists is known a priori. To specify the banded alignment
 * the developer can use the seqan3::align_cfg::band_fixed_size option.
 * This band configuration is initialised with a seqan3::align_cfg::lower_diagonal and an
 * seqan3::align_cfg::upper_diagonal. The upper diagonal must always be greater than or equal to the lower diagonal.
 * To choose the correct band parameters, imagine a matrix with the first sequence written on top and the second sequence
 * along the left vertical side. A negative value reflects a start of the diagonal within the vertical part while a
 * positive value implies a start within the top part of this matrix at the respective position.
 *
 * ## Global and local alignments
 *
 * The standard global and local alignments can be configured using seqan3::align_cfg::method_global and .
 * seqan3::align_cfg::method_local, respectively.
 * The global alignment can be further refined using the seqan3::align_cfg::sequence_ends configuration.
 * This configuration allows to enable, respectively disable the scoring of leading and trailing gaps at the respective
 * sequence ends. This option is not available for the local alignment where scoring gaps at the ends of the sequences
 * is always disabled.
 *
 * ## Algorithmic details
 *
 * By default a generic alignment algorithm is used that supports all valid alignment configurations, but for some
 * special combinations of parameters a notably faster algorithm is available.
 * It is automatically selected if all of the following requirements are satisfied:
 *  * Edit distance gaps, i.e. seqan3::align_cfg::gap is initialised with default initialised seqan3::gap_scheme
 *  * Edit distance scoring for \ref nucleotide "nucleotide alphabets", i.e. seqan3::align_cfg::scoring is initialised with default initialised seqan3::nucleotide_scoring_scheme.
 *  * Global alignment, i.e. seqan3::align_cfg::method_global.
 *
 * There is a special shortcut for the above required scoring/gap configs called seqan3::align_cfg::edit_scheme,
 * which can be used to safe some typing.
 *
 * The edit configuration can be further specialised with following configs:
 *  * Allow maximal number of errors, i.e specify the seqan3::align_cfg::max_error configuration
 *  * Compute a semi-global alignment, i.e seqan3::align_cfg::aligned_ends is initialised with seqan3::end_gaps::free_ends_first.
 *
 * \note If there was a configuration that is not suitable for the edit distance algorithm the standard alignment
 * algorithm is executed as a fallback.
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
