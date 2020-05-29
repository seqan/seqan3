# The SeqAn Cookbook {#cookbook}

[TOC]

This document provides example recipes on how to carry out particular tasks using the SeqAn functionalities in C++.
Please note that these **recipes are not ordered**. You can use the links in the table of contents or the search
function of your browser to navigate them.

It will take some time, but we hope to expand this document into containing numerous great examples.
If you have suggestions for how to improve the Cookbook and/or examples you would like included,
please feel free to contact us.

#  Read sequence files
\snippet doc/tutorial/introduction/introduction_file_input.cpp fileinput

# Write a custom validator
This recipe implements a validator that checks if a numeric argument is an integral square (i.e.  0, 1, 4, 9...).
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
  3. prints the string and it's reverse complement
  4. prints the six-frame translation of the string

\include doc/tutorial/ranges/range_solution3.cpp

# Reading records
After construction, you can now read the sequence records.
Our file object behaves like a range so you can use a range based for loop to conveniently iterate over the file:

\include test/snippet/io/sequence_file/sequence_file_input_record_iter.cpp

\attention An input file is a **single input range**, which means you can only iterate over it **once**!

\note It is important to write `auto &` and not just `auto`, otherwise you will copy the record on every iteration.

You can also use structured binding, i.e. `for (auto & [seq, id, qual] : fin)`  
**But beware: with structured bindings you do need to get the order of elements correctly!**

You can also read a file in chunks:

# Reading records in chunks

\snippet doc/tutorial/sequence_file/sequence_file_snippets.cpp read_in_batches

The example above will iterate over the file by reading 10 records at a time.
If no 10 records are available any more, it will just print the remaining records

# Applying a filter to a file
In some occasions you are only interested in sequence records that fulfill a certain criterion,
e.g. having a minimum sequence length or a minimum average quality.

This recipe can be used to filter the sequences in your file on a custom criterion.

\snippet doc/tutorial/sequence_file/sequence_file_snippets.cpp quality_filter

# Reading paired-end reads
In modern Next Generation Sequencing experiments you often have paired-end read data which is split into two files.
The read pairs are identified by their identical name/id and that their position in the two files is the same.

This recipe can be used to handle one pair of reads at a time.

\snippet doc/tutorial/sequence_file/sequence_file_snippets.cpp paired_reads

# Storing records in a std::vector
This recipe creates a small program that reads in a FASTA file and stores all the records in a std::vector.

\snippet doc/tutorial/sequence_file/sequence_file_solution2.cpp solution

Note that you can move the record out of the file if you want to store it somewhere without copying.

\snippet doc/tutorial/sequence_file/sequence_file_snippets.cpp record_type2

# Writing records

The easiest way to write to a sequence file is to use the seqan3::sequence_file_output::push_back()
or seqan3::sequence_file_output::emplace_back() member functions.
These work similarly to how they work on an std::vector.

\include test/snippet/io/sequence_file/sequence_file_output_record_wise_iteration.cpp

# File conversion

\snippet doc/tutorial/sequence_file/sequence_file_snippets.cpp file_conversion

# Define a custom scoring scheme

\snippet doc/tutorial/pairwise_alignment/configurations.cpp include_scoring_scheme

\snippet doc/tutorial/pairwise_alignment/configurations.cpp scoring_scheme

\attention SeqAn's alignment algorithm computes the maximal similarity score, thus the match score must be set to a
positive value and the score for mismatch and gaps must be negative in order to maximize over the matching letters.

# Calculate edit distance for a set of sequences

This recipe can be used to calculate the edit distance for all six pairwise combinations. Here we only allow at most 7
errors and filter all alignments with 6 or less errors.

\include doc/tutorial/pairwise_alignment/pairwise_alignment_solution_6.cpp

# Searching for matches

This recipe can be used to search for all occurrences of a substring and print the number of hits and the
positions in an ascending ordering.

\include doc/tutorial/search/search_solution2.cpp

If you want to allow for errors in your query, you need to configure the approximate search with a
seqan3::search_cfg::max_error or seqan3::search_cfg::max_error_rate object.

To search for either 1 insertion or 1 deletion you can use the seqan3::search_cfg::max_error:

\snippet doc/tutorial/search/search_small_snippets.cpp error_search

# Reading the CIGAR information into an actual alignment

In SeqAn, the conversion from a CIGAR string to an alignment (two aligned_sequences) is done automatically for you.
You can access it by querying seqan3::field::alignment from the record:

\snippet doc/tutorial/alignment_file/alignment_file_snippets.cpp alignments_without_ref

# Combining sequence and alignment files

This recipe can be used to:
1. Read in a FASTA file with the reference and a SAM file with the alignment
2. Filter the alignment records and only take those with a mapping quality >= 30.
3. For the resulting alignments, print which read was mapped against the reference id and the number of seqan3::gap's
   involved in the alignment (either in aligned reference or in read sequence).

\snippet doc/tutorial/alignment_file/alignment_file_solution2.cpp solution

# Map reads ans write output to SAM file
For a full recipe on creating your own readmapper, see the very end of the tutorial \ref tutorial_read_mapper.

\snippet doc/tutorial/read_mapper/read_mapper_step4.cpp solution

# Constructing a basic argument parser

\include doc/tutorial/read_mapper/read_mapper_indexer_step1.cpp

# Constructing a subcommand argument parser

\include doc/howto/subcommand_argument_parser/subcommand_arg_parse.cpp

#  Serialise data structures with cereal

\include doc/howto/use_cereal/load.hpp
