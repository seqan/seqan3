# The SeqAn Cookbook {#cookbook}

<!-- SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
     SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
     SPDX-License-Identifier: CC-BY-4.0
-->

[TOC]

This document provides example recipes on how to carry out particular tasks using the SeqAn functionalities in C++.
Please note that these **recipes are not ordered**. You can use the links in the table of contents or the search
function of your browser to navigate them.

It will take some time, but we hope to expand this document into containing numerous great examples.
If you have suggestions for how to improve the Cookbook and/or examples you would like included,
please feel free to contact us.

# Read sequence files
\snippet doc/cookbook/file_input.cpp fileinput

# Construction and assignment of alphabet symbols

\snippet doc/tutorial/04_alphabet/alphabet_main.cpp create
\snippet doc/tutorial/04_alphabet/alphabet_main.cpp closing
\snippet doc/tutorial/04_alphabet/alphabet_main.cpp rank

# Reverse complement and the six-frame translation of a string using views
This recipe creates a small program that
  1. reads a string from the command line (first argument to the program)
  2. "converts" the string to a range of seqan3::dna5 (Bonus: throws an exception if loss of information occurs)
  3. prints the string and its reverse complement
  4. prints the six-frame translation of the string

\include doc/tutorial/05_ranges/range_solution3.cpp

# Reading records
After construction, you can now read the sequence records.
Our file object behaves like a range, you can use a range-based for loop to conveniently iterate over the file:

\include test/snippet/io/sequence_file/sequence_file_input_record_iter.cpp

\attention An input file is a **single input range**, which means you can only iterate over it **once**!

\note It is important to write `auto &` and not just `auto`, otherwise you will copy the record on every iteration.

You can also use structured binding, i.e. `for (auto & [seq, id, qual] : fin)`
**But beware: with structured bindings you do need to get the order of elements correct!**

You can also read a file in chunks:

# Reading records in chunks

\snippet doc/tutorial/07_sequence_file/sequence_file_read_in_batches.cpp main

The example above will iterate over the file by reading 10 records at a time.
If no 10 records are available anymore, it will just print the remaining records.

# Applying a filter to a file
On some occasions you are only interested in sequence records that fulfill a certain criterion,
e.g. having a minimum sequence length or a minimum average quality.

This recipe can be used to filter the sequences in your file by a custom criterion.

\snippet doc/tutorial/07_sequence_file/sequence_file_quality_filter.cpp main

# Reading paired-end reads
In modern Next Generation Sequencing experiments you often have paired-end read data which is split into two files.
The read pairs are identified by their identical name/id and position in the two files.

This recipe can be used to handle one pair of reads at a time.

\snippet doc/tutorial/07_sequence_file/sequence_file_paired_reads.cpp main

# Storing records in a std::vector
This recipe creates a small program that reads in a FASTA file and stores all the records in a std::vector.

\snippet doc/tutorial/07_sequence_file/sequence_file_solution2.cpp solution

Note that you can move the record out of the file if you want to store it somewhere without copying.

\snippet doc/tutorial/07_sequence_file/sequence_file_move_record.cpp main

# Writing records

The easiest way to write to a sequence file is to use the seqan3::sequence_file_output::push_back()
or seqan3::sequence_file_output::emplace_back() member functions.
These work similarly to how they work on a std::vector.

\include doc/tutorial/07_sequence_file/sequence_file_output_record.cpp

\cond DEV
The class seqan3::sequence_file_output takes an extra parameter allowing to custom select the
fields and their order.

\include test/snippet/io/sequence_file/sequence_file_output_fields_trait_1.cpp
\endcond

# File conversion

\snippet doc/tutorial/07_sequence_file/sequence_file_file_conversion.cpp main

# Define a custom scoring scheme

\snippet doc/tutorial/08_pairwise_alignment/configurations.cpp include_scoring_scheme

\snippet doc/tutorial/08_pairwise_alignment/configurations.cpp scoring_scheme

\attention SeqAn's alignment algorithm computes the maximal similarity score, thus the match score must be set to a
positive value and the scores for mismatch and gap must be negative in order to maximize over the matching letters.

# Calculate edit distance for a set of sequences

This recipe can be used to calculate the edit distance for all six pairwise combinations. Here we only allow at most 7
errors and filter all alignments with 6 or fewer errors.

\include doc/tutorial/08_pairwise_alignment/pairwise_alignment_solution_6.cpp

# Searching for matches

This recipe can be used to search for all occurrences of a substring and print the number of hits and the
positions in an ascending ordering.

\include doc/tutorial/09_search/search_solution2.cpp

If you want to allow errors in your query, you need to configure the approximate search with the following search
configuration objects:
- seqan3::search_cfg::max_error_total: Maximum number of total errors
- seqan3::search_cfg::max_error_substitution: Maximum number of substitutions
- seqan3::search_cfg::max_error_insertion: Maximum number of insertions
- seqan3::search_cfg::max_error_deletion: Maximum number of deletions
These are constructed with absolute numbers or rates:
- seqan3::search_cfg::error_count: Absolute number of errors
- seqan3::search_cfg::error_rate: Rate of errors \f$\in[0,1]\f$

To search for either 1 insertion or 1 deletion you can use the seqan3::search_cfg::error_count:

\snippet doc/tutorial/09_search/search_small_snippets.cpp error_search

# Reading the CIGAR information from a SAM file and constructing an alignment

This recipe can be used to:
1. Read in a FASTA file with the reference and a SAM file with the alignment
2. Filter the alignment records and only take those with a mapping quality >= 30.
3. For the resulting alignments, print which read was mapped against with reference id and the number of seqan3::gap's
   involved in the alignment (either in aligned reference or in read sequence).

\snippet doc/tutorial/10_sam_file/sam_file_solution2.cpp solution

# Map reads and write output to SAM file
For a full recipe on creating your own readmapper, see the very end of the tutorial \ref tutorial_read_mapper.

\snippet doc/tutorial/11_read_mapper/read_mapper_step4.cpp solution

# Constructing a basic argument parser

\include doc/tutorial/11_read_mapper/read_mapper_indexer_step1.cpp

# Constructing a subcommand argument parser

\include doc/howto/subcommand_argument_parser/subcommand_arg_parse.cpp

# Serialise data structures with cereal

\include doc/howto/use_cereal/load.hpp

# Converting a range of an alphabet {#cookbook_convert_alphabet_range}

\include doc/cookbook/alphabet_conversion.cpp

# A custom dna4 alphabet that converts all unknown characters to A {#cookbook_custom_dna4_alphabet}

When assigning from `char` or converting from a larger nucleotide alphabet to a smaller one, *loss of information*
can occur since obviously some bases are not available. When converting to seqan3::dna5 or seqan3::rna5,
non-canonical bases (letters other than A, C, G, T, U) are converted to ``'N'`` to preserve ambiguity at that position.
For seqan3::dna4 and seqan3::rna4 there is no letter ``'N'`` to represent ambiguity, so the conversion from `char` for
IUPAC characters tries to choose the best fitting alternative (see seqan3::dna4 for more details).

If you would like to always convert unknown characters to `A` instead, you can create your own alphabet with a
respective char conversion table very easily like this:

\include doc/cookbook/custom_dna4.cpp

If you are interested in custom alphabets, also take a look at our tutorial \ref howto_write_an_alphabet.

# Controlling threads of (de-)compression streams {#setting_compression_threads}

When \ref io_compression "reading or writing compressed files", parallelisation is automatically applied when
using BGZF-compressed files, e.g., BAM files.
This will use `4` threads by default and can be adjusted by setting `seqan3::contrib::bgzf_thread_count` to
the desired value:

\snippet doc/cookbook/compression_threads.cpp example

# Auto vectorized dna4 complement

Our alphabet seqan3::dna4 cannot be easily auto-vectorized by the compiler.

See [this discussion](https://github.com/seqan/seqan3/issues/1970) for more details.

You can add your own alphabet that is auto-vectorizable in some use cases.
Here is an example for a dna4-like alphabet:

\snippet test/performance/simd_dna4.hpp cookbook

# All SeqAn documentation snippets

The following lists all snippets that appear in our documentation.
Search for keywords with `Strg + F`.

<!-- ALL SNIPPETS START -->
\include test/snippet/alignment/cigar_conversion/alignment_from_cigar.cpp
\include test/snippet/alignment/cigar_conversion/alignment_from_cigar_io.cpp
\include test/snippet/alignment/cigar_conversion/cigar_from_alignment.cpp
\include test/snippet/alignment/cigar_conversion/cigar_from_alignment_with_clipping.cpp
\include test/snippet/alignment/configuration/align_cfg_band_example.cpp
\include test/snippet/alignment/configuration/align_cfg_edit_example.cpp
\include test/snippet/alignment/configuration/align_cfg_gap_cost_affine_example.cpp
\include test/snippet/alignment/configuration/align_cfg_method_global.cpp
\include test/snippet/alignment/configuration/align_cfg_method_local.cpp
\include test/snippet/alignment/configuration/align_cfg_min_score_example.cpp
\include test/snippet/alignment/configuration/align_cfg_on_result.cpp
\include test/snippet/alignment/configuration/align_cfg_output_alignment.cpp
\include test/snippet/alignment/configuration/align_cfg_output_begin_position.cpp
\include test/snippet/alignment/configuration/align_cfg_output_end_position.cpp
\include test/snippet/alignment/configuration/align_cfg_output_examples.cpp
\include test/snippet/alignment/configuration/align_cfg_output_score.cpp
\include test/snippet/alignment/configuration/align_cfg_output_sequence1_id.cpp
\include test/snippet/alignment/configuration/align_cfg_output_sequence2_id.cpp
\include test/snippet/alignment/configuration/align_cfg_parallel_example.cpp
\include test/snippet/alignment/configuration/align_cfg_score_type.cpp
\include test/snippet/alignment/configuration/align_cfg_vectorised_example.cpp
\include test/snippet/alignment/configuration/minimal_alignment_config.cpp
\include test/snippet/alignment/matrix/detail/alignment_matrix_column_major_range_base.cpp
\include test/snippet/alignment/matrix/detail/debug_matrix_score.cpp
\include test/snippet/alignment/matrix/detail/debug_matrix_trace.cpp
\include test/snippet/alignment/pairwise/align_pairwise.cpp
\include test/snippet/alignment/pairwise/align_pairwise_range.cpp
\include test/snippet/alignment/pairwise/alignment_configurator.cpp
\include test/snippet/alignment/pairwise/parallel_align_pairwise_with_callback.cpp
\include test/snippet/alignment/scoring/aminoacid_scoring_scheme.cpp
\include test/snippet/alignment/scoring/nucleotide_scoring_scheme.cpp
\include test/snippet/alphabet/all_ambiguous.cpp
\include test/snippet/alphabet/all_literal.cpp
\include test/snippet/alphabet/all_nonambiguous.cpp
\include test/snippet/alphabet/alphabet_base.cpp
\include test/snippet/alphabet/alphabet_size.cpp
\include test/snippet/alphabet/aminoacid/aa10li.cpp
\include test/snippet/alphabet/aminoacid/aa10li_char_literal.cpp
\include test/snippet/alphabet/aminoacid/aa10li_literal.cpp
\include test/snippet/alphabet/aminoacid/aa10murphy.cpp
\include test/snippet/alphabet/aminoacid/aa10murphy_char_literal.cpp
\include test/snippet/alphabet/aminoacid/aa10murphy_literal.cpp
\include test/snippet/alphabet/aminoacid/aa20.cpp
\include test/snippet/alphabet/aminoacid/aa20_char_literal.cpp
\include test/snippet/alphabet/aminoacid/aa20_literal.cpp
\include test/snippet/alphabet/aminoacid/aa27.cpp
\include test/snippet/alphabet/aminoacid/aa27_char_literal.cpp
\include test/snippet/alphabet/aminoacid/aa27_literal.cpp
\include test/snippet/alphabet/aminoacid/enable_aminoacid.cpp
\include test/snippet/alphabet/assign_char_strictly_to.cpp
\include test/snippet/alphabet/assign_char_to.cpp
\include test/snippet/alphabet/assign_rank_to.cpp
\include test/snippet/alphabet/char_is_valid_for.cpp
\include test/snippet/alphabet/cigar/cigar.cpp
\include test/snippet/alphabet/cigar/cigar_assign_string.cpp
\include test/snippet/alphabet/cigar/cigar_get_index.cpp
\include test/snippet/alphabet/cigar/cigar_get_type.cpp
\include test/snippet/alphabet/cigar/cigar_operation.cpp
\include test/snippet/alphabet/cigar/cigar_operation_char_literal.cpp
\include test/snippet/alphabet/cigar/cigar_value_assignment.cpp
\include test/snippet/alphabet/cigar/cigar_value_construction.cpp
\include test/snippet/alphabet/composite/alphabet_tuple_base_subtype_assignment.cpp
\include test/snippet/alphabet/composite/alphabet_tuple_base_subtype_construction.cpp
\include test/snippet/alphabet/composite/alphabet_tuple_base_value_assignment.cpp
\include test/snippet/alphabet/composite/alphabet_tuple_base_value_construction.cpp
\include test/snippet/alphabet/composite/alphabet_variant.cpp
\include test/snippet/alphabet/composite/alphabet_variant_char_representation.cpp
\include test/snippet/alphabet/composite/alphabet_variant_conversion.cpp
\include test/snippet/alphabet/composite/alphabet_variant_conversion_explicit.cpp
\include test/snippet/alphabet/composite/alphabet_variant_is_alternative.cpp
\include test/snippet/alphabet/composite/alphabet_variant_subtype_construction.cpp
\include test/snippet/alphabet/composite/alphabet_variant_value_construction.cpp
\include test/snippet/alphabet/composite/semialphabet_any.cpp
\include test/snippet/alphabet/container/bitpacked_sequence.cpp
\include test/snippet/alphabet/container/concatenated_sequences.cpp
\include test/snippet/alphabet/container/concatenated_sequences_insert.cpp
\include test/snippet/alphabet/container/concatenated_sequences_insert2.cpp
\include test/snippet/alphabet/gap/gap.cpp
\include test/snippet/alphabet/gap/gapped.cpp
\include test/snippet/alphabet/mask/mask.cpp
\include test/snippet/alphabet/mask/masked.cpp
\include test/snippet/alphabet/nucleotide/complement_cpo.cpp
\include test/snippet/alphabet/nucleotide/dna15.cpp
\include test/snippet/alphabet/nucleotide/dna15_char_literal.cpp
\include test/snippet/alphabet/nucleotide/dna15_implicit_conversion_from_rna15.cpp
\include test/snippet/alphabet/nucleotide/dna15_implicit_conversion_from_rna15_inherit.cpp
\include test/snippet/alphabet/nucleotide/dna15_implicit_conversion_from_rna15_vector.cpp
\include test/snippet/alphabet/nucleotide/dna15_implicit_conversion_from_rna15_views.cpp
\include test/snippet/alphabet/nucleotide/dna15_literal.cpp
\include test/snippet/alphabet/nucleotide/dna16sam.cpp
\include test/snippet/alphabet/nucleotide/dna16sam_char_literal.cpp
\include test/snippet/alphabet/nucleotide/dna16sam_literal.cpp
\include test/snippet/alphabet/nucleotide/dna3bs.cpp
\include test/snippet/alphabet/nucleotide/dna3bs_char_literal.cpp
\include test/snippet/alphabet/nucleotide/dna3bs_literal.cpp
\include test/snippet/alphabet/nucleotide/dna4.cpp
\include test/snippet/alphabet/nucleotide/dna4_char_literal.cpp
\include test/snippet/alphabet/nucleotide/dna4_implicit_conversion_from_rna4.cpp
\include test/snippet/alphabet/nucleotide/dna4_implicit_conversion_from_rna4_inherit.cpp
\include test/snippet/alphabet/nucleotide/dna4_implicit_conversion_from_rna4_vector.cpp
\include test/snippet/alphabet/nucleotide/dna4_implicit_conversion_from_rna4_views.cpp
\include test/snippet/alphabet/nucleotide/dna4_literal.cpp
\include test/snippet/alphabet/nucleotide/dna5.cpp
\include test/snippet/alphabet/nucleotide/dna5_char_literal.cpp
\include test/snippet/alphabet/nucleotide/dna5_implicit_conversion_from_rna5.cpp
\include test/snippet/alphabet/nucleotide/dna5_implicit_conversion_from_rna5_inherit.cpp
\include test/snippet/alphabet/nucleotide/dna5_implicit_conversion_from_rna5_vector.cpp
\include test/snippet/alphabet/nucleotide/dna5_implicit_conversion_from_rna5_views.cpp
\include test/snippet/alphabet/nucleotide/dna5_literal.cpp
\include test/snippet/alphabet/nucleotide/rna15.cpp
\include test/snippet/alphabet/nucleotide/rna15_char_literal.cpp
\include test/snippet/alphabet/nucleotide/rna15_implicit_conversion_from_dna15.cpp
\include test/snippet/alphabet/nucleotide/rna15_implicit_conversion_from_dna15_inherit.cpp
\include test/snippet/alphabet/nucleotide/rna15_implicit_conversion_from_dna15_vector.cpp
\include test/snippet/alphabet/nucleotide/rna15_implicit_conversion_from_dna15_views.cpp
\include test/snippet/alphabet/nucleotide/rna15_literal.cpp
\include test/snippet/alphabet/nucleotide/rna4.cpp
\include test/snippet/alphabet/nucleotide/rna4_char_literal.cpp
\include test/snippet/alphabet/nucleotide/rna4_implicit_conversion_from_dna4.cpp
\include test/snippet/alphabet/nucleotide/rna4_implicit_conversion_from_dna4_inherit.cpp
\include test/snippet/alphabet/nucleotide/rna4_implicit_conversion_from_dna4_vector.cpp
\include test/snippet/alphabet/nucleotide/rna4_implicit_conversion_from_dna4_views.cpp
\include test/snippet/alphabet/nucleotide/rna4_literal.cpp
\include test/snippet/alphabet/nucleotide/rna5.cpp
\include test/snippet/alphabet/nucleotide/rna5_char_literal.cpp
\include test/snippet/alphabet/nucleotide/rna5_implicit_conversion_from_dna5.cpp
\include test/snippet/alphabet/nucleotide/rna5_implicit_conversion_from_dna5_inherit.cpp
\include test/snippet/alphabet/nucleotide/rna5_implicit_conversion_from_dna5_vector.cpp
\include test/snippet/alphabet/nucleotide/rna5_implicit_conversion_from_dna5_views.cpp
\include test/snippet/alphabet/nucleotide/rna5_literal.cpp
\include test/snippet/alphabet/quality/phred42.cpp
\include test/snippet/alphabet/quality/phred42_char_literal.cpp
\include test/snippet/alphabet/quality/phred42_literal.cpp
\include test/snippet/alphabet/quality/phred63.cpp
\include test/snippet/alphabet/quality/phred63_char_literal.cpp
\include test/snippet/alphabet/quality/phred63_literal.cpp
\include test/snippet/alphabet/quality/phred68solexa.cpp
\include test/snippet/alphabet/quality/phred68solexa_char_literal.cpp
\include test/snippet/alphabet/quality/phred68solexa_literal.cpp
\include test/snippet/alphabet/quality/phred94.cpp
\include test/snippet/alphabet/quality/phred94_char_literal.cpp
\include test/snippet/alphabet/quality/phred94_literal.cpp
\include test/snippet/alphabet/quality/qualified.cpp
\include test/snippet/alphabet/structure/dot_bracket3.cpp
\include test/snippet/alphabet/structure/dot_bracket3_char_literal.cpp
\include test/snippet/alphabet/structure/dot_bracket3_literal.cpp
\include test/snippet/alphabet/structure/dssp9.cpp
\include test/snippet/alphabet/structure/dssp9_char_literal.cpp
\include test/snippet/alphabet/structure/dssp9_literal.cpp
\include test/snippet/alphabet/structure/structured_aa.cpp
\include test/snippet/alphabet/structure/structured_rna.cpp
\include test/snippet/alphabet/structure/wuss.cpp
\include test/snippet/alphabet/structure/wuss_char_literal.cpp
\include test/snippet/alphabet/structure/wuss_is_pair_close.cpp
\include test/snippet/alphabet/structure/wuss_is_pair_open.cpp
\include test/snippet/alphabet/structure/wuss_is_unpaired.cpp
\include test/snippet/alphabet/structure/wuss_literal.cpp
\include test/snippet/alphabet/structure/wuss_max_pseudoknot_depth.cpp
\include test/snippet/alphabet/structure/wuss_pseudoknot_id.cpp
\include test/snippet/alphabet/views/char_strictly_to.cpp
\include test/snippet/alphabet/views/char_to.cpp
\include test/snippet/alphabet/views/complement.cpp
\include test/snippet/alphabet/views/range_view_to_char.cpp
\include test/snippet/alphabet/views/range_view_to_rank.cpp
\include test/snippet/alphabet/views/rank_to.cpp
\include test/snippet/alphabet/views/to_char.cpp
\include test/snippet/alphabet/views/to_rank.cpp
\include test/snippet/alphabet/views/translate_dna5.cpp
\include test/snippet/alphabet/views/translate_join.cpp
\include test/snippet/alphabet/views/translate_usage.cpp
\include test/snippet/alphabet/views/trim_quality_dna5q.cpp
\include test/snippet/alphabet/views/trim_quality_phred42.cpp
\include test/snippet/alphabet/views/validate_char_for.cpp
\include test/snippet/argument_parser/argument_parser_1.cpp
\include test/snippet/argument_parser/argument_parser_2.cpp
\include test/snippet/argument_parser/argument_parser_3.cpp
\include test/snippet/argument_parser/auxiliary.cpp
\include test/snippet/argument_parser/custom_argument_parsing_enumeration.cpp
\include test/snippet/argument_parser/custom_enumeration.cpp
\include test/snippet/argument_parser/is_option_set.cpp
\include test/snippet/argument_parser/validators_1.cpp
\include test/snippet/argument_parser/validators_2.cpp
\include test/snippet/argument_parser/validators_3.cpp
\include test/snippet/argument_parser/validators_4.cpp
\include test/snippet/argument_parser/validators_chaining.cpp
\include test/snippet/argument_parser/validators_input_directory.cpp
\include test/snippet/argument_parser/validators_input_file.cpp
\include test/snippet/argument_parser/validators_input_file_ext_from_file.cpp
\include test/snippet/argument_parser/validators_output_directory.cpp
\include test/snippet/argument_parser/validators_output_file.cpp
\include test/snippet/argument_parser/validators_output_file_ext_from_file.cpp
\include test/snippet/core/add_enum_bitwise_operators.cpp
\include test/snippet/core/cereal_example.cpp
\include test/snippet/core/configuration/configuration_combine.cpp
\include test/snippet/core/configuration/configuration_compatibility.cpp
\include test/snippet/core/configuration/configuration_get.cpp
\include test/snippet/core/configuration/configuration_get_or.cpp
\include test/snippet/core/configuration/configuration_setup.cpp
\include test/snippet/core/debug_stream_flags.cpp
\include test/snippet/core/debug_stream_set_underlying_stream.cpp
\include test/snippet/core/debug_stream_set_underlying_stream2.cpp
\include test/snippet/core/debug_stream_usage.cpp
\include test/snippet/core/detail/customisation_point.cpp
\include test/snippet/core/detail/deferred_crtp_base.cpp
\include test/snippet/core/detail/is_class_template_declarable_with.cpp
\include test/snippet/core/detail/strong_type_adding_skills.cpp
\include test/snippet/core/detail/strong_type_error_window.cpp
\include test/snippet/core/detail/strong_type_new_usage.cpp
\include test/snippet/core/detail/strong_type_usage.cpp
\include test/snippet/core/detail/template_inspection_usage.cpp
\include test/snippet/core/detail/template_inspection_usage_2.cpp
\include test/snippet/core/detail/template_inspection_usage_3.cpp
\include test/snippet/io/detail/detail_record.cpp
\include test/snippet/io/detail/iterator_write_range.cpp
\include test/snippet/io/detail/safe_filesystem_entry_snippet.cpp
\include test/snippet/io/record_1.cpp
\include test/snippet/io/record_2.cpp
\include test/snippet/io/sam_file/begin_iterator.cpp
\include test/snippet/io/sam_file/emplace_back.cpp
\include test/snippet/io/sam_file/push_back_record.cpp
\include test/snippet/io/sam_file/push_back_tuple.cpp
\include test/snippet/io/sam_file/record_based_writing.cpp
\include test/snippet/io/sam_file/record_based_writing2.cpp
\include test/snippet/io/sam_file/sam_file_input_begin_and_front.cpp
\include test/snippet/io/sam_file/sam_file_input_construction_from_filename.cpp
\include test/snippet/io/sam_file/sam_file_input_construction_from_stream.cpp
\include test/snippet/io/sam_file/sam_file_input_construction_without_automatic_type_deduction.cpp
\include test/snippet/io/sam_file/sam_file_input_front.cpp
\include test/snippet/io/sam_file/sam_file_input_get_header.cpp
\include test/snippet/io/sam_file/sam_file_input_my_traits.cpp
\include test/snippet/io/sam_file/sam_file_input_reading_custom_fields.cpp
\include test/snippet/io/sam_file/sam_file_input_reading_filter.cpp
\include test/snippet/io/sam_file/sam_file_input_reading_move_record.cpp
\include test/snippet/io/sam_file/sam_file_input_reading_range_based_for_loop.cpp
\include test/snippet/io/sam_file/sam_file_input_reading_structured_bindings.cpp
\include test/snippet/io/sam_file/sam_file_output_cout_write.cpp
\include test/snippet/io/sam_file/sam_file_output_custom_fields.cpp
\include test/snippet/io/sam_file/sam_file_output_filename_construction.cpp
\include test/snippet/io/sam_file/sam_file_output_filename_construction_with_ref_info.cpp
\include test/snippet/io/sam_file/sam_file_output_format_construction.cpp
\include test/snippet/io/sam_file/sam_file_output_input_range.cpp
\include test/snippet/io/sam_file/sam_file_output_io_pipeline.cpp
\include test/snippet/io/sam_file/sam_file_output_set_header.cpp
\include test/snippet/io/sam_file/sam_file_output_write_range.cpp
\include test/snippet/io/sam_file/sam_flags.cpp
\include test/snippet/io/sam_file/sam_tag_dictionary/general_usage.cpp
\include test/snippet/io/sam_file/sam_tag_dictionary/sam_tag_dictionary.cpp
\include test/snippet/io/sam_file/sam_tag_dictionary/unknown_tag.cpp
\include test/snippet/io/sequence_file/sequence_file_input_aminoacid.cpp
\include test/snippet/io/sequence_file/sequence_file_input_auto_ref.cpp
\include test/snippet/io/sequence_file/sequence_file_input_custom_fields.cpp
\include test/snippet/io/sequence_file/sequence_file_input_decomposed.cpp
\include test/snippet/io/sequence_file/sequence_file_input_file_view.cpp
\include test/snippet/io/sequence_file/sequence_file_input_istringstream.cpp
\include test/snippet/io/sequence_file/sequence_file_input_record_iter.cpp
\include test/snippet/io/sequence_file/sequence_file_input_record_move.cpp
\include test/snippet/io/sequence_file/sequence_file_input_return_record.cpp
\include test/snippet/io/sequence_file/sequence_file_input_template_deduction.cpp
\include test/snippet/io/sequence_file/sequence_file_input_template_specification.cpp
\include test/snippet/io/sequence_file/sequence_file_input_trait_overwrite.cpp
\include test/snippet/io/sequence_file/sequence_file_output_batch_write.cpp
\include test/snippet/io/sequence_file/sequence_file_output_col_based_writing.cpp
\include test/snippet/io/sequence_file/sequence_file_output_cout_write.cpp
\include test/snippet/io/sequence_file/sequence_file_output_direct_writing.cpp
\include test/snippet/io/sequence_file/sequence_file_output_emplace_back.cpp
\include test/snippet/io/sequence_file/sequence_file_output_fields_trait_1.cpp
\include test/snippet/io/sequence_file/sequence_file_output_fields_trait_2.cpp
\include test/snippet/io/sequence_file/sequence_file_output_push_back_record.cpp
\include test/snippet/io/sequence_file/sequence_file_output_push_back_tuple.cpp
\include test/snippet/io/sequence_file/sequence_file_output_range_interface.cpp
\include test/snippet/io/sequence_file/sequence_file_output_record_wise_iteration.cpp
\include test/snippet/io/sequence_file/sequence_file_output_template_deduction.cpp
\include test/snippet/io/sequence_file/sequence_file_output_view_pipeline.cpp
\include test/snippet/io/structure_file/structure_file_input_arg_spec.cpp
\include test/snippet/io/structure_file/structure_file_input_auto_temp_deduc.cpp
\include test/snippet/io/structure_file/structure_file_input_data_out.cpp
\include test/snippet/io/structure_file/structure_file_input_filter_criteria.cpp
\include test/snippet/io/structure_file/structure_file_input_mod_traits.cpp
\include test/snippet/io/structure_file/structure_file_input_move.cpp
\include test/snippet/io/structure_file/structure_file_input_record_iter.cpp
\include test/snippet/io/structure_file/structure_file_input_ref_return.cpp
\include test/snippet/io/structure_file/structure_file_input_skip_fields.cpp
\include test/snippet/io/structure_file/structure_file_input_structured_bindings.cpp
\include test/snippet/io/structure_file/structure_file_input_trait_def.cpp
\include test/snippet/io/structure_file/structure_file_output_col_based.cpp
\include test/snippet/io/structure_file/structure_file_output_emplace_back.cpp
\include test/snippet/io/structure_file/structure_file_output_equal.cpp
\include test/snippet/io/structure_file/structure_file_output_iter_by_rec.cpp
\include test/snippet/io/structure_file/structure_file_output_mult_rec.cpp
\include test/snippet/io/structure_file/structure_file_output_pass_rec.cpp
\include test/snippet/io/structure_file/structure_file_output_pipe_func.cpp
\include test/snippet/io/structure_file/structure_file_output_pipeline.cpp
\include test/snippet/io/structure_file/structure_file_output_push_back.cpp
\include test/snippet/io/structure_file/structure_file_output_push_back_2.cpp
\include test/snippet/io/structure_file/structure_file_output_temp_param_deduc.cpp
\include test/snippet/io/structure_file/structure_file_output_write_fields.cpp
\include test/snippet/io/structure_file/structure_file_output_write_std_out.cpp
\include test/snippet/io/views/async_input_buffer.cpp
\include test/snippet/io/views/detail/take_exactly.cpp
\include test/snippet/io/views/detail/take_line_view_adaptor_def.cpp
\include test/snippet/io/views/detail/take_line_view_behaviour.cpp
\include test/snippet/io/views/detail/take_line_view_tokenise.cpp
\include test/snippet/io/views/detail/take_until_view.cpp
\include test/snippet/range/detail/random_access_iterator.cpp
\include test/snippet/range/views/range_view_all_composability.cpp
\include test/snippet/range/views/range_view_all_notation.cpp
\include test/snippet/range/views/range_view_all_retransform.cpp
\include test/snippet/release/3.0.0_complement-and-translate-without-copying.cpp
\include test/snippet/release/3.0.0_file-conversion.cpp
\include test/snippet/release/3.0.0_read-fasta-and-print-entities.cpp
\include test/snippet/release/3.0.1_parallel-alignments.cpp
\include test/snippet/release/3.0.2_dynamic-hit-configuration.cpp
\include test/snippet/release/3.0.2_lazy-search-result-range.cpp
\include test/snippet/release/3.0.2_user-callback-in-alignment-and-search.cpp
\include test/snippet/release/3.0.3_io-record.cpp
\include test/snippet/release/3.0.3_literals-namespace.cpp
\include test/snippet/search/bi_fm_index.cpp
\include test/snippet/search/bi_fm_index_collection.cpp
\include test/snippet/search/bi_fm_index_cursor_cycle.cpp
\include test/snippet/search/bi_fm_index_cursor_extend_left_seq.cpp
\include test/snippet/search/bi_fm_index_cursor_to_fwd_cursor.cpp
\include test/snippet/search/configuration_default.cpp
\include test/snippet/search/configuration_error.cpp
\include test/snippet/search/configuration_on_result.cpp
\include test/snippet/search/configuration_output.cpp
\include test/snippet/search/configuration_parallel.cpp
\include test/snippet/search/dream_index/counting_agent.cpp
\include test/snippet/search/dream_index/counting_agent_construction.cpp
\include test/snippet/search/dream_index/counting_vector.cpp
\include test/snippet/search/dream_index/example_query_genome_region.cpp
\include test/snippet/search/dream_index/interleaved_bloom_filter_clear.cpp
\include test/snippet/search/dream_index/interleaved_bloom_filter_constructor.cpp
\include test/snippet/search/dream_index/interleaved_bloom_filter_constructor_compressed.cpp
\include test/snippet/search/dream_index/interleaved_bloom_filter_constructor_uncompress.cpp
\include test/snippet/search/dream_index/interleaved_bloom_filter_emplace.cpp
\include test/snippet/search/dream_index/interleaved_bloom_filter_increase_bin_number_to.cpp
\include test/snippet/search/dream_index/membership_agent_bulk_contains.cpp
\include test/snippet/search/dream_index/membership_agent_construction.cpp
\include test/snippet/search/dynamic_hit_configuration_example.cpp
\include test/snippet/search/fm_index.cpp
\include test/snippet/search/fm_index_collection.cpp
\include test/snippet/search/fm_index_cursor.cpp
\include test/snippet/search/hit_configuration_examples.cpp
\include test/snippet/search/kmer_index/shape.cpp
\include test/snippet/search/search.cpp
\include test/snippet/search/search_with_user_callback.cpp
\include test/snippet/search/views/kmer_hash.cpp
\include test/snippet/search/views/minimiser.cpp
\include test/snippet/search/views/minimiser_hash.cpp
\include test/snippet/snippet_main.cpp
\include test/snippet/std/view/subrange.cpp
\include test/snippet/test/simd_utility.cpp
\include test/snippet/test/tmp_directory.cpp
\include test/snippet/utility/bloom_filter/bloom_filter_constructor.cpp
\include test/snippet/utility/bloom_filter/bloom_filter_constructor_compressed.cpp
\include test/snippet/utility/bloom_filter/bloom_filter_contains.cpp
\include test/snippet/utility/bloom_filter/bloom_filter_count.cpp
\include test/snippet/utility/bloom_filter/bloom_filter_emplace.cpp
\include test/snippet/utility/bloom_filter/bloom_filter_reset.cpp
\include test/snippet/utility/char_operations/char_predicate.cpp
\include test/snippet/utility/char_operations/char_predicate_is_char.cpp
\include test/snippet/utility/char_operations/char_predicate_is_in_interval.cpp
\include test/snippet/utility/char_operations/char_predicate_is_operator.cpp
\include test/snippet/utility/container/aligned_allocator.cpp
\include test/snippet/utility/container/dynamic_bitset_begin.cpp
\include test/snippet/utility/container/dynamic_bitset_binary_and_member.cpp
\include test/snippet/utility/container/dynamic_bitset_binary_or_member.cpp
\include test/snippet/utility/container/dynamic_bitset_binary_xor_member.cpp
\include test/snippet/utility/container/dynamic_bitset_construct_string.cpp
\include test/snippet/utility/container/dynamic_bitset_construct_uint64_t.cpp
\include test/snippet/utility/container/dynamic_bitset_flip.cpp
\include test/snippet/utility/container/dynamic_bitset_flip_all.cpp
\include test/snippet/utility/container/dynamic_bitset_flip_pos.cpp
\include test/snippet/utility/container/dynamic_bitset_reset_all.cpp
\include test/snippet/utility/container/dynamic_bitset_reset_pos.cpp
\include test/snippet/utility/container/dynamic_bitset_resize.cpp
\include test/snippet/utility/container/dynamic_bitset_set_all.cpp
\include test/snippet/utility/container/dynamic_bitset_set_pos.cpp
\include test/snippet/utility/container/dynamic_bitset_shift_left_copy.cpp
\include test/snippet/utility/container/dynamic_bitset_shift_left_inplace.cpp
\include test/snippet/utility/container/dynamic_bitset_shift_right_copy.cpp
\include test/snippet/utility/container/dynamic_bitset_shift_right_inplace.cpp
\include test/snippet/utility/container/dynamic_bitset_subscript.cpp
\include test/snippet/utility/container/dynamic_bitset_usage.cpp
\include test/snippet/utility/container/small_string.cpp
\include test/snippet/utility/detail/ceil_log2.cpp
\include test/snippet/utility/detail/floor_log2.cpp
\include test/snippet/utility/pow.cpp
\include test/snippet/utility/simd/detail/builtin_simd.cpp
\include test/snippet/utility/simd/detail/default_simd_max_length.cpp
\include test/snippet/utility/simd/detail/is_builtin_simd.cpp
\include test/snippet/utility/simd/fill.cpp
\include test/snippet/utility/simd/iota.cpp
\include test/snippet/utility/simd/simd.cpp
\include test/snippet/utility/simd/simd_extract.cpp
\include test/snippet/utility/simd/simd_load.cpp
\include test/snippet/utility/simd/simd_store.cpp
\include test/snippet/utility/simd/simd_traits.cpp
\include test/snippet/utility/simd/simd_transpose.cpp
\include test/snippet/utility/simd/simd_upcast.cpp
\include test/snippet/utility/simd/views/iota_simd.cpp
\include test/snippet/utility/simd/views/iota_simd_transform.cpp
\include test/snippet/utility/simd/views/to_simd.cpp
\include test/snippet/utility/tuple/pod_tuple.cpp
\include test/snippet/utility/tuple_utility.cpp
\include test/snippet/utility/type_list/detail/type_list_algorithm_all_of.cpp
\include test/snippet/utility/type_list/detail/type_list_algorithm_for_each.cpp
\include test/snippet/utility/type_list/list_traits_at.cpp
\include test/snippet/utility/type_list/list_traits_back.cpp
\include test/snippet/utility/type_list/list_traits_concat.cpp
\include test/snippet/utility/type_list/list_traits_contains.cpp
\include test/snippet/utility/type_list/list_traits_count.cpp
\include test/snippet/utility/type_list/list_traits_drop.cpp
\include test/snippet/utility/type_list/list_traits_drop_front.cpp
\include test/snippet/utility/type_list/list_traits_drop_last.cpp
\include test/snippet/utility/type_list/list_traits_find.cpp
\include test/snippet/utility/type_list/list_traits_find_if.cpp
\include test/snippet/utility/type_list/list_traits_front.cpp
\include test/snippet/utility/type_list/list_traits_replace_at.cpp
\include test/snippet/utility/type_list/list_traits_size.cpp
\include test/snippet/utility/type_list/list_traits_split_after.cpp
\include test/snippet/utility/type_list/list_traits_take.cpp
\include test/snippet/utility/type_list/list_traits_take_last.cpp
\include test/snippet/utility/type_list/list_traits_transform.cpp
\include test/snippet/utility/type_pack/detail/type_pack_algorithm_all_of.cpp
\include test/snippet/utility/type_pack/detail/type_pack_algorithm_for_each.cpp
\include test/snippet/utility/type_pack/pack_traits_at.cpp
\include test/snippet/utility/type_pack/pack_traits_back.cpp
\include test/snippet/utility/type_pack/pack_traits_count.cpp
\include test/snippet/utility/type_pack/pack_traits_drop.cpp
\include test/snippet/utility/type_pack/pack_traits_drop_front.cpp
\include test/snippet/utility/type_pack/pack_traits_drop_last.cpp
\include test/snippet/utility/type_pack/pack_traits_find.cpp
\include test/snippet/utility/type_pack/pack_traits_find_if.cpp
\include test/snippet/utility/type_pack/pack_traits_front.cpp
\include test/snippet/utility/type_pack/pack_traits_replace_at.cpp
\include test/snippet/utility/type_pack/pack_traits_size.cpp
\include test/snippet/utility/type_pack/pack_traits_split_after.cpp
\include test/snippet/utility/type_pack/pack_traits_take.cpp
\include test/snippet/utility/type_pack/pack_traits_take_last.cpp
\include test/snippet/utility/type_pack/pack_traits_transform.cpp
\include test/snippet/utility/type_traits/function_traits.cpp
\include test/snippet/utility/type_traits/lazy_conditional.cpp
\include test/snippet/utility/type_traits/transformation_trait_or.cpp
\include test/snippet/utility/views/chunk.cpp
\include test/snippet/utility/views/convert_15_to_5.cpp
\include test/snippet/utility/views/convert_int_to_bool.cpp
\include test/snippet/utility/views/deep_no_param.cpp
\include test/snippet/utility/views/deep_pass_ref.cpp
\include test/snippet/utility/views/deep_with_param.cpp
\include test/snippet/utility/views/deep_wrap_args.cpp
\include test/snippet/utility/views/elements.cpp
\include test/snippet/utility/views/enforce_random_access.cpp
\include test/snippet/utility/views/interleave.cpp
\include test/snippet/utility/views/pairwise_combine.cpp
\include test/snippet/utility/views/repeat.cpp
\include test/snippet/utility/views/repeat_n.cpp
\include test/snippet/utility/views/single_pass_input.cpp
\include test/snippet/utility/views/slice.cpp
\include test/snippet/utility/views/type_reduce.cpp
