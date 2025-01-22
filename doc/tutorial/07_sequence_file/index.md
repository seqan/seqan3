# Sequence File Input and Output {#tutorial_sequence_file}

<!-- SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
     SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
     SPDX-License-Identifier: CC-BY-4.0
-->

<b>Learning Objective:</b> <br/>
You will get an overview of how Input/Output files are handled in SeqAn and learn how to read and write
sequence files. This tutorial is a walk-through with links into the API documentation and also meant as a
source for copy-and-paste code.

\tutorial_head{Easy, 90 min, \ref setup "Setup"\, \ref tutorial_ranges\, \ref tutorial_alphabets "Alphabets",}

[TOC]

# File I/O in SeqAn

Most file formats in bioinformatics are structured as lists of records.
In SeqAn we model our files as a **range over records**.
This interface allows us to easily stream over a file, apply filters and convert formats,
sometimes only in a single line of code.
The file format is automatically detected by the file name extension and compressed files can be handled
without any effort. We can even stream over files in a python-like way with range-based for loops:

\snippet doc/tutorial/07_sequence_file/sequence_file_record_range.cpp record_range

We will explain the details about reading and writing files in the [Sequence File](#section_sequence_files) section below.
Currently, SeqAn supports the following file formats:

- Sequence file formats:
  - seqan3::format_fasta
  - seqan3::format_fastq
  - seqan3::format_embl
  - seqan3::format_sam (alignment format contains enough information to be read/written as pure sequence file)

- Structure file formats:
  - seqan3::format_vienna

- Alignment file formats:
  - seqan3::format_sam

\warning Access to compressed files relies on external libraries.
For instance, you need to have *zlib* installed for reading `.gz` files and *libbz2* for reading `.bz2` files.
You can check whether you have installed these libraries by running `cmake .` in your build directory.
If `-- Optional dependency: ZLIB-x.x.x found.` is displayed on the command line then you can read/write
compressed files in your programs.

## Basic layout of SeqAn file objects

Before we dive into the details, we will outline the general design of our file objects,
hoping that it will make the following tutorial easier to understand.

As mentioned above, our file object is a range over records. More specifically, a range over objects of type
seqan3::sequence_record.

Output files can handle various types that fulfill the requirements of the format (e.g.
a sequence has to be a range over an alphabet).
In contrast to this, input files have certain default types for record fields that can be modified via a *traits type*.

For example, on construction you can specify seqan3::sequence_file_input_default_traits_dna or
seqan3::sequence_file_input_default_traits_aa for reading `dna` (\ref alphabet_nucleotide) and `protein` (\ref alphabet_aminoacid)
sequences respectively (section [File traits](#section_file_traits) will covers this in more detail).

Opening and closing files is also handled automatically.
If a file cannot be opened for reading or writing, a seqan3::file_open_error is thrown.

# Sequence file formats

Sequence files are the most generic and common biological files.
Well-known formats include FASTA and FASTQ.

### FASTA format

A FASTA record contains the sequence id and the sequence characters. Here is an example of a FASTA file:

\include doc/tutorial/07_sequence_file/example.fasta

In SeqAn we provide the seqan3::format_fasta to read sequence files in FASTA format.

### FASTQ format

A FASTQ record contains an additional quality value for each sequence character. Here is an example of a FASTQ file:

\include doc/tutorial/07_sequence_file/example.fastq

In SeqAn we provide the seqan3::format_fastq to read sequence files in FASTQ format.

### EMBL format

An EMBL record stores sequence and its annotation together. We are only interested in the id (ID) and sequence (SQ)
information for a sequence file. Qualities are not stored in this format.
Here is an example of an EMBL file:

\include doc/tutorial/07_sequence_file/example.embl

In SeqAn we provide the seqan3::format_embl to read sequence files in EMBL format.

### File extensions

The formerly introduced formats can be identified by the following file name extensions
(this is important for automatic format detection from a file name as you will learn in the next section).

| File Format | SeqAn format class   | File Extensions                                   |
| ------------| ---------------------|---------------------------------------------------|
| FASTA       | seqan3::format_fasta |   `.fa`, `.fasta`, `.fna`, `.ffn`, `.ffa`, `.frn` |
| FASTQ       | seqan3::format_fastq |   `.fq`, `.fastq`                                 |
| EMBL        | seqan3::format_embl  |   `.embl`                                         |


You can access the valid file extensions via the `file_extensions` member variable in a format
and you can also customise this list if you want to allow different or additional file extensions:

\snippet doc/tutorial/07_sequence_file/sequence_file_file_extensions.cpp main

# Fields {#section_sequence_files}

The Sequence file abstraction supports reading three different fields:

  1. seqan3::sequence_record::sequence
  2. seqan3::sequence_record::id
  3. seqan3::sequence_record::base_qualities

The three fields are retrieved by default.

# Reading a sequence file

You can include the SeqAn sequence file functionality with:

\snippet doc/tutorial/07_sequence_file/sequence_file_record_range.cpp include

## Construction

At first, you need to construct a seqan3::sequence_file_input object that handles file access.
In most cases, you construct from a file name, but you can also construct a sequence file object directly from a stream
(e.g. std::cin or std::stringstream), but then you need to know your format beforehand:

\snippet doc/tutorial/07_sequence_file/sequence_file_filename_construction.cpp main

All template parameters of the seqan3::sequence_file_input are automatically deduced, even the format!
We **detect the format by the file name extension**.
The file extension in the example above is `.fasta` so we choose the seqan3::format_fasta.


### File traits {#section_file_traits}

The seqan3::sequence_file_input needs to know the types of the data you are reading in
(e.g. that your sequence is a `dna` sequence) at **compile-time**.
These necessary types are defined in the  **traits type** argument (seqan3::sequence_file_input::traits_type)
which by default is set to seqan3::sequence_file_input_default_traits_dna.

We thereby assume that
* you want to read the sequence into a std::vector over a seqan3::dna5 alphabet,
* store the sequence names/id in a std::string and
* the qualities in a std::vector over the seqan3::phred42 alphabet.

In case you want to read a **protein** sequence instead, we also provide the
seqan3::sequence_file_input_default_traits_aa traits type which sets the seqan3::sequence_record::sequence field to a
std::vector over the seqan3::aa27 alphabet.

You can specify the traits object as the first template argument for the sequence file:

\snippet doc/tutorial/07_sequence_file/sequence_file_amino_acid_type_trait.cpp main

You can also customise the types by inheriting from one of the default traits and changing the type manually.
See the detailed information on seqan3::sequence_file_input_default_traits_dna for an example.

## Reading records

After construction, you can now read the seqan3::sequence_record's. As described in the basic file layout,
our file objects behave like ranges so you can use a range based for loop to conveniently iterate over the file:

\include test/snippet/io/sequence_file/sequence_file_input_record_iter.cpp

\attention An input file is a **single input range**, which means you can only iterate over it **once**!

In the above example, `record` has the type seqan3::sequence_file_input::record_type which is seqan3::sequence_record.

\note It is important to write `auto &` and not just `auto`, otherwise you will copy the record on every iteration.

You can read up more on the different ways to stream over the file object in the detailed documentation
of seqan3::sequence_file_input.

\assignment{Assignment 1: Reading a FASTQ file}
Copy and paste the following FASTQ file to some location, e.g. the tmp directory:

\snippet doc/tutorial/07_sequence_file/sequence_file_solution1.cpp fastq_file

and then create a simple program that reads in all the records from that file and prints the id,
sequence and quality information to the command line using the seqan3::debug_stream.
Access the record via the member accessor of seqan3::sequence_record!

Note: you include the seqan3::debug_stream with the following header

\snippet doc/tutorial/07_sequence_file/sequence_file_solution1.cpp include_debug_stream

\endassignment
\solution

\snippet doc/tutorial/07_sequence_file/sequence_file_solution1.cpp solution

The code will print the following:
\include doc/tutorial/07_sequence_file/sequence_file_solution1.err
\endsolution

## The record type

In the examples above, we always use `auto` to deduce the record type automatically.
In case you need the type explicitly, e.g. if you want to store the records in a variable,
you can access the type member seqan3::sequence_file_input::record_type.

\snippet doc/tutorial/07_sequence_file/sequence_file_record_type.cpp main

You can move the record out of the file if you want to store it somewhere without copying.

\snippet doc/tutorial/07_sequence_file/sequence_file_move_record.cpp main

\assignment{Assignment 2: Storing records in a std::vector}

Create a small program that reads in a FASTA file and stores all the records in a std::vector.

After reading, print the vector (this works natively with the seqan3::debug_stream).

Test your program with the following file:

\snippet doc/tutorial/07_sequence_file/sequence_file_solution2.cpp fasta_file

It should print the following:

\include doc/tutorial/07_sequence_file/sequence_file_solution2.err
Note that the quality (third element in the triple) is empty because we are reading a FASTA file.
\endassignment
\solution
\snippet doc/tutorial/07_sequence_file/sequence_file_solution2.cpp solution
\endsolution

# Sequence files as views

Since SeqAn files are ranges, you can also create views over files.
This enables us to create solutions for lots of use cases using only a few lines of code.

## Reading a file in chunks

A common use case is to read chunks from a file instead of the whole file at once or line by line.

You can do so easily on a file range by using the seqan3::views::chunk.

\snippet doc/tutorial/07_sequence_file/sequence_file_read_in_batches.cpp main

The example above will iterate over the file by reading 10 records at a time.
If no 10 records are available any more, it will just print the remaining records.

## Applying a filter to a file

On some occasions, you are only interested in sequence records that fulfill a certain criteria,
e.g. having a minimum sequence length or a minimum average quality.
Just like in the example with *seqan3::views::chunk* you can use *std::ranges::filter* for this purpose:

\snippet doc/tutorial/07_sequence_file/sequence_file_quality_filter.cpp main

To remind you what you have learned in the \ref tutorial_ranges tutorial before,
a view is not applied immediately but **lazy evaluated**.
That means your file is still parsed record by record and not at once.

## Reading paired-end reads {#sequence_file_section_fun_with_ranges}

In modern Next Generation Sequencing experiments, you often encounter paired-end read data which is split into two files.
The read pairs are identified by their identical id and position in the files.

If you want to handle one pair of reads at a time, you can do so easily with a seqan3::views::zip.

\snippet doc/tutorial/07_sequence_file/sequence_file_paired_reads.cpp main

\assignment{Assignment 3: Fun with file ranges}

Implement a small program that reads in a FASTQ file and prints the first 2 sequences
that have a length of at least 5.

Hints:
* You can use `std::ranges::size` to retrieve the size of a range.
* You need the following includes for `std::views::filter` and `std::views::take`
\snippet doc/tutorial/07_sequence_file/sequence_file_solution3.cpp include_ranges


Test your program on the following FASTQ file:

\snippet doc/tutorial/07_sequence_file/sequence_file_solution3.cpp fastq_file

It should output `seq1` and `seq3`.

\endassignment
\solution
\snippet doc/tutorial/07_sequence_file/sequence_file_solution3.cpp solution
\endsolution

# Writing a sequence file

## Construction

You construct the seqan3::sequence_file_output just like the seqan3::sequence_file_input
by giving it a file name or a stream.

\snippet test/snippet/io/sequence_file/sequence_file_output_template_deduction.cpp main

Writing to std::cout:

\include test/snippet/io/sequence_file/sequence_file_output_cout_write.cpp

## Writing records

The easiest way to write to a sequence file is to use the seqan3::sequence_file_output::push_back() member function.
It works similarly to how it works on a std::vector.

\include doc/tutorial/07_sequence_file/sequence_file_output_record.cpp

\assignment{Assignment 4: Writing a FASTQ file}

Use your code (or the solution) from the previous assignment.
Iterate over the records with a for loop and instead of just printing the ids,
write out **all** the records that satisfy the filter to a new file called `output.fastq`.

Test your code on the same FASTQ file.
The file `output.fastq` should contain the following records:

\snippet doc/tutorial/07_sequence_file/sequence_file_solution4.cpp fastq_file

\endassignment
\solution
\snippet doc/tutorial/07_sequence_file/sequence_file_solution4.cpp solution
\endsolution

## Files as views

Again we want to point out the convenient advantage of modelling files as ranges. The seqan3::sequence_file_input models
std::ranges::input_range and in the "reading a file" section you already saw a few examples of how to pipe a view onto
such a range. The output files model std::ranges::output_range. This allows us to use ranges algorithms like
std::ranges::move to "move" records from input to output files.

\cond DEV
\todo There is no standard way for piping output ranges (e.g. `fin | fout;`), yet. Reevaluate this section later.

In the same way you can pipe the output file:
\endcond

\snippet doc/tutorial/07_sequence_file/sequence_file_piping_in_out.cpp main

\assignment{Assignment 5: Fun with file ranges 2}

Working on your solution from the previous assignment, try to remove the for loop in favour of a std::ranges algorithm.

The result should be the same.

\endassignment
\solution
\snippet doc/tutorial/07_sequence_file/sequence_file_solution5.cpp solution
\endsolution

# File conversion

As mentioned before, the seqan3::sequence_record is shared between formats which allows for easy file conversion.
For example you can read in a FASTQ file and output a FASTA file in one line:

\snippet doc/tutorial/07_sequence_file/sequence_file_file_conversion.cpp main

Yes, that's it! Of course this only works because all fields that are required in FASTA are provided in FASTQ.
The other way around would not work as easily because we have no quality information (and would make less sense too).
