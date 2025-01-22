<!-- SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
     SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
     SPDX-License-Identifier: CC-BY-4.0
-->

### Writing Sequence Files

Sequence files are the most generic and common biological files. Well-known formats include
FASTA and FASTQ, but some may also be interested in treating SAM or BAM files as sequence
files, discarding the alignment.

The Sequence file abstraction supports writing three different fields:
  1. seqan3::field::seq
  2. seqan3::field::id
  3. seqan3::field::qual

The member functions take any and either of these fields. If the field ID of an argument cannot be deduced, it
is assumed to correspond to the field ID of the respective template parameter.

#### Construction and specialisation

This class comes with two constructors, one for construction from a file name and one for construction from
an existing stream and a known format. The first one automatically picks the format based on the extension
of the file name. The second can be used if you have a non-file stream, like std::cout or std::ostringstream,
that you want to read from and/or if you cannot use file-extension based detection, but know that your output
file has a certain format.
<br><br>
In most cases the template parameters are deduced completely automatically:

\snippet test/snippet/io/sequence_file/sequence_file_output_template_deduction.cpp main

Writing to std::cout:

\include test/snippet/io/sequence_file/sequence_file_output_cout_write.cpp

Note that this is not the same as writing `sequence_file_output<>` (with angle brackets). In the latter case they are
explicitly set to their default values, in the former case
[automatic deduction](https://en.cppreference.com/w/cpp/language/class_template_argument_deduction) happens which
chooses different parameters depending on the constructor arguments. Prefer deduction over explicit defaults.

#### Writing record-wise

You can iterate over this file record-wise:

\include test/snippet/io/sequence_file/sequence_file_output_record_wise_iteration.cpp

The easiest way to write to a sequence file is to use the push_back() or emplace_back() member functions. These
work similarly to how they work on a std::vector. If you pass a tuple to push_back() or give arguments to
emplace_back() the seqan3::field ID of the i-th tuple-element/argument is assumed to be the i-th value of
selected_field_ids, i.e. by default the first is assumed to be seqan3::field::seq, the second seqan3::field::id
and the third one seqan3::field::qual. You may give less fields than are selected if the actual format you are
writing to can cope with less
(e.g. for FASTA it is sufficient to write seqan3::field::seq and seqan3::field::id, even if selected_field_ids
also contains seqan3::field::qual at the third position).
You may also use the output file's iterator for writing, however, this rarely provides an advantage.

#### Writing record-wise (custom fields)

If you want to change the order of the parameters, you can pass a non-empty fields trait object to the
sequence_file_output constructor to select the fields that are used for interpreting the arguments.
The following snippets demonstrates the usage of such a fields trait object.

\include test/snippet/io/sequence_file/sequence_file_output_fields_trait_1.cpp

A different way of passing custom fields to the file is to pass a seqan3::record – instead of a tuple – to
push_back(). The seqan3::record clearly indicates which of its elements has which seqan3::field ID so the file will
use that information instead of the template argument. This is especially handy when reading from one file and
writing to another, because you don't have to configure the output file to match the input file, it will just work:

\include test/snippet/io/sequence_file/sequence_file_output_fields_trait_2.cpp

#### Writing record-wise in batches

You can write multiple records at once, by assigning to the file:

\include test/snippet/io/sequence_file/sequence_file_output_batch_write.cpp

#### File I/O pipelines

Record-wise writing in batches also works for writing from input files directly to output files, because input
files are also input ranges in SeqAn:

\include test/snippet/io/sequence_file/sequence_file_output_direct_writing.cpp

This can be combined with file-based views to create I/O pipelines:

\snippet test/snippet/io/sequence_file/sequence_file_output_view_pipeline.cpp snippet

#### Column-based writing

The record-based interface treats the file as a range of tuples (the records), but in certain situations
you might have the data as columns, i.e. a tuple-of-ranges, instead of range-of-tuples.
You can use column-based writing in that case, it uses operator=() and seqan3::views::zip():

\include test/snippet/io/sequence_file/sequence_file_output_col_based_writing.cpp

#### Formats

We currently support writing the following formats:
	* seqan3::format_fasta
	* seqan3::format_fastq
	* seqan3::format_embl
	* seqan3::format_genbank
	* seqan3::format_sam
