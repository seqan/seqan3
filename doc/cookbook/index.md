# The SeqAn Cookbook {#cookbook}

[TOC]

This document provides example recipes on how to carry out particular tasks using the SeqAn functionalities in C++.
Please note that these **recipes are not ordered**. You can use the links in the table of contents or the search
function of your browser to navigate them.

It will take some time, but we hope to expand this document into containing numerous great examples.
If you have suggestions for how to improve the Cookbook and/or examples you would like included,
please feel free to contact us.

# Read sequence files
\snippet doc/tutorial/introduction/introduction_file_input.cpp fileinput

# Write a custom validator
This recipe implements a validator that checks if a numeric argument is an integral square (i.e. 0, 1, 4, 9...).
Invalid values throw a seqan3::validation_error.

\snippet doc/tutorial/concepts/custom_validator_solution2.cpp validator

# Construction and assignment of alphabet symbols

\snippet doc/tutorial/alphabet/alphabet_main.cpp create
\snippet doc/tutorial/alphabet/alphabet_main.cpp closing
\snippet doc/tutorial/alphabet/alphabet_main.cpp rank

# Reverse complement and the six-frame translation of a string using views
This recipe creates a small program that
  1. reads a string from the command line (first argument to the program)
  2. "converts" the string to a range of seqan3::dna5 (Bonus: throws an exception if loss of information occurs)
  3. prints the string and its reverse complement
  4. prints the six-frame translation of the string

\include doc/tutorial/ranges/range_solution3.cpp

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

\snippet doc/tutorial/sequence_file/sequence_file_read_in_batches.cpp main

The example above will iterate over the file by reading 10 records at a time.
If no 10 records are available anymore, it will just print the remaining records.

# Applying a filter to a file
On some occasions you are only interested in sequence records that fulfill a certain criterion,
e.g. having a minimum sequence length or a minimum average quality.

This recipe can be used to filter the sequences in your file by a custom criterion.

\snippet doc/tutorial/sequence_file/sequence_file_quality_filter.cpp main

# Reading paired-end reads
In modern Next Generation Sequencing experiments you often have paired-end read data which is split into two files.
The read pairs are identified by their identical name/id and position in the two files.

This recipe can be used to handle one pair of reads at a time.

\snippet doc/tutorial/sequence_file/sequence_file_paired_reads.cpp main

# Storing records in a std::vector
This recipe creates a small program that reads in a FASTA file and stores all the records in a std::vector.

\snippet doc/tutorial/sequence_file/sequence_file_solution2.cpp solution

Note that you can move the record out of the file if you want to store it somewhere without copying.

\snippet doc/tutorial/sequence_file/sequence_file_move_record.cpp main

# Writing records

The easiest way to write to a sequence file is to use the seqan3::sequence_file_output::push_back()
or seqan3::sequence_file_output::emplace_back() member functions.
These work similarly to how they work on a std::vector.

\include doc/tutorial/sequence_file//sequence_file_output_record.cpp

\cond DEV
The class seqan3::sequence_file_output takes an extra parameter allowing to custom select the
fields and their order.

\include test/snippet/io/sequence_file/sequence_file_output_fields_trait_1.cpp
\endcond

# File conversion

\snippet doc/tutorial/sequence_file/sequence_file_file_conversion.cpp main

# Define a custom scoring scheme

\snippet doc/tutorial/pairwise_alignment/configurations.cpp include_scoring_scheme

\snippet doc/tutorial/pairwise_alignment/configurations.cpp scoring_scheme

\attention SeqAn's alignment algorithm computes the maximal similarity score, thus the match score must be set to a
positive value and the scores for mismatch and gap must be negative in order to maximize over the matching letters.

# Calculate edit distance for a set of sequences

This recipe can be used to calculate the edit distance for all six pairwise combinations. Here we only allow at most 7
errors and filter all alignments with 6 or fewer errors.

\include doc/tutorial/pairwise_alignment/pairwise_alignment_solution_6.cpp

# Searching for matches

This recipe can be used to search for all occurrences of a substring and print the number of hits and the
positions in an ascending ordering.

\include doc/tutorial/search/search_solution2.cpp

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

\snippet doc/tutorial/search/search_small_snippets.cpp error_search

# Reading the CIGAR information into an actual alignment

In SeqAn, the conversion from a CIGAR string to an alignment (two aligned_sequences) is done automatically for you.
You can access it by querying seqan3::field::alignment from the record:

\snippet doc/tutorial/sam_file/sam_file_snippets.cpp alignments_without_ref

# Combining sequence and alignment files

This recipe can be used to:
1. Read in a FASTA file with the reference and a SAM file with the alignment
2. Filter the alignment records and only take those with a mapping quality >= 30.
3. For the resulting alignments, print which read was mapped against with reference id and the number of seqan3::gap's
   involved in the alignment (either in aligned reference or in read sequence).

\snippet doc/tutorial/sam_file/sam_file_solution2.cpp solution

# Map reads and write output to SAM file
For a full recipe on creating your own readmapper, see the very end of the tutorial \ref tutorial_read_mapper.

\snippet doc/tutorial/read_mapper/read_mapper_step4.cpp solution

# Constructing a basic argument parser

\include doc/tutorial/read_mapper/read_mapper_indexer_step1.cpp

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
This will use all available threads by default and can be adjusted by setting `seqan3::contrib::bgzf_thread_count` to
the desired value:

\snippet doc/cookbook/compression_threads.cpp example
