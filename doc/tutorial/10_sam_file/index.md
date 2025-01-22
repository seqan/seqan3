# SAM Input and Output in SeqAn {#tutorial_sam_file}

<!-- SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
     SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
     SPDX-License-Identifier: CC-BY-4.0
-->

***Learning Objective:***

<b>Learning Objective:</b> <br/>
You will get an overview of how to read and write SAM/BAM files.
This tutorial is a walk-through with links into the API documentation and is also meant as a source for copy-and-paste code.

\tutorial_head{Medium, 60 min, \ref setup\, \ref tutorial_alphabets\, \ref tutorial_sequence_file,}

[TOC]

# Introduction

SAM files are used to store pairwise alignments between two (biological) sequences. There are also other output formats,
like BLAST, which can store sequence alignments, but in this tutorial we will focus on SAM/BAM files.
In addition to the alignment, these formats store information such as the start positions or mapping qualities.
SAM files are a little more complex than sequence files but the basic design is the same.
If you are new to SeqAn, we strongly recommend completing the tutorial \ref tutorial_sequence_file first.

# SAM/BAM file formats

## SAM format

SAM stands for Sequence Alignment/Map format. It is a TAB-delimited text format consisting of a header
section, which is optional, and an alignment section
(see the [official SAM specifications](https://samtools.github.io/hts-specs/SAMv1.pdf)).

Here is an example of a SAM file:

\include doc/tutorial/10_sam_file/example.sam

The following table summarises the columns of a SAM file:

| #  | SAM Column ID |  Description                                      |
|:--:|:--------------|:--------------------------------------------------|
| 1  | QNAME         | Query template NAME                               |
| 2  | FLAG          | bitwise FLAG                                      |
| 3  | RNAME         | Reference sequence NAME                           |
| 4  | POS           | 1-based leftmost mapping POSition                 |
| 5  | MAPQ          | MAPping Quality                                   |
| 6  | CIGAR         | CIGAR string                                      |
| 7  | RNEXT         | Reference name of the mate/next read              |
| 8  | PNEXT         | Position of the mate/next read                    |
| 9  | TLEN          | observed Template LENgth                          |
| 10 | SEQ           | segment SEQuence                                  |
| 11 | QUAL          | ASCII of Phred-scaled base QUALity+33             |

If you want to read more about the SAM format, take a look at the
[official specifications](https://samtools.github.io/hts-specs/SAMv1.pdf).

## BAM format

BAM is the binary format version of SAM. It provides the same data as the SAM format with negligible and subtle
differences in most use cases.

# SAM file fields

To make things clearer, here is the table of SAM columns and the corresponding fields of a SAM file record:

| #  | SAM Column ID |  FIELD name                            | seqan3::field             |
|:--:|:--------------|:---------------------------------------|:------------------------- |
| 1  | QNAME         | seqan3::sam_record::id                 | seqan3::field::id         |
| 2  | FLAG          | seqan3::sam_record::flag               | seqan3::field::flag       |
| 3  | RNAME         | seqan3::sam_record::reference_id       | seqan3::field::ref_id     |
| 4  | POS           | seqan3::sam_record::reference_position | seqan3::field::ref_offset |
| 5  | MAPQ          | seqan3::sam_record::mapping_quality    | seqan3::field::mapq       |
| 6  | CIGAR         | seqan3::sam_record::cigar_sequence     | seqan3::field::cigar      |
| 7  | RNEXT         | seqan3::sam_record::mate_reference_id  | seqan3::field::mate       |
| 8  | PNEXT         | seqan3::sam_record::mate_position      | seqan3::field::mate       |
| 9  | TLEN          | seqan3::sam_record::template_length    | seqan3::field::mate       |
| 10 | SEQ           | seqan3::sam_record::sequence           | seqan3::field::seq        |
| 11 | QUAL          | seqan3::sam_record::base_qualities     | seqan3::field::qual       |

SAM files provide following additional fields:
* seqan3::sam_record::tags (seqan3::field::tags)
* seqan3::sam_record::header_ptr (seqan3::field::header_ptr)

## File extensions

The formerly introduced formats can be identified by the following file name extensions
(this is important for automatic format detection from a file name as you will learn in the next section).

| File Format  | File Extensions   |
| -------------|-------------------|
| SAM          |   `.sam`          |
| BAM          |   `.bam`          |

You can access and modify the valid file extensions via the `file_extension` member variable in a format tag:

\snippet doc/tutorial/10_sam_file/sam_file_file_extensions.cpp main

# Reading SAM files

Before we start, you should copy and paste this [example file](example.sam) into a file location of your choice
(we use the current path in the examples, so make sure you adjust your path).

\attention Make sure the file you copied is tab delimited!

## Construction

The construction works analogously to sequence files by passing a file name,
in which case, all template parameters are automatically deduced (by the file name extension).
Or you can pass a stream (e.g. std::cin or std::stringstream), but then you need to know your format beforehand:

\snippet doc/tutorial/10_sam_file/sam_file_filename_construction.cpp main

## Accessing individual record members

You can access a record member like this:

\snippet doc/tutorial/10_sam_file/sam_file_sam_record.cpp main

See seqan3::sam_record for all data accessors.

\assignment{Assignment 1: Accumulating mapping qualities}

Let's assume we want to compute the average mapping quality of a SAM file.

For this purpose, write a small program that
    * only reads the mapping quality (seqan3::sam_record::mapping_quality) of a SAM file and
    * computes the average of all qualities.

Use the following file to test your program:

\snippet doc/tutorial/10_sam_file/sam_file_solution1.cpp sam_file

It should output:
```
Average: 27.4
```

\endassignment
\solution

\snippet doc/tutorial/10_sam_file/sam_file_solution1.cpp solution

\endsolution

# Alignment representation in SAM/BAM files

The SAM format is the common output format of read mappers where you align short read sequences
to one or more large reference sequences.
In fact, the SAM format stores those alignment information only partially:
It **does not store the reference sequence** but only the query/read sequence and a *CIGAR* string
representing the alignment based on the read.

Take this SAM record as an example:

```
r003   73   ref   3   17   1M1D4M   *  0   0  TAGGC   *
```

The record gives you the following information:
A read with name `r003` has been mapped to a reference with name `ref` at position `3`
(in the reference, counting from 1) with a quality of `17` (Phred scaled).
The flag has a value of `73` which indicates that the read is paired, the first in pair, but the mate is unmapped
(see [this website](https://broadinstitute.github.io/picard/explain-flags.html) for a nice explanation of SAM flags).
Fields set to `0` or `*` indicate empty fields and contain no valuable information.

The cigar string is `1M1D4M` which represents the following alignment:

```
      1 2 3 4 5 6 7 8 9 ...
ref   N N N N N N N N N ...
read      T - A G G C
```

where the reference sequence is not known (represented by `N`).
You will learn in the next section how to handle additional reference sequence information.

If you want to read up more about cigar strings,
take a look at the [SAM specifications](https://samtools.github.io/hts-specs/SAMv1.pdf)
or the [SAMtools paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2723002/).

## Reading the CIGAR string

By default, the seqan3::sam_file_input will always read the seqan3::sam_record::cigar_sequence and store it
into a std::vector\<seqan3::cigar\>:

\snippet doc/tutorial/10_sam_file/sam_file_read_cigar.cpp main

## Transforming the CIGAR information into an alignment

In SeqAn, we represent an alignment as a tuple of two `seqan3::aligned_sequence`s.

The conversion from a CIGAR string to an alignment can be done with the function `seqan3::alignment_from_cigar`.
You need to pass the reference sequence with the position the read was aligned to and the read sequence:

\include snippet/alignment/cigar_conversion/alignment_from_cigar_io.cpp

The code will print the following:
\include snippet/alignment/cigar_conversion/alignment_from_cigar_io.err

\assignment{Assignment 2: Combining sequence and alignment files}

Read the following reference sequence FASTA file (see the sequence file tutorial if you need a reminder):

\snippet doc/tutorial/10_sam_file/sam_file_solution2.cpp ref_file

Then read the following SAM file while providing the reference sequence information.

\snippet doc/tutorial/10_sam_file/sam_file_solution2.cpp sam_file

Only use
* seqan3::sam_record::id,
* seqan3::sam_record::reference_id,
* seqan3::sam_record::mapping_quality, and
* seqan3::sam_record::cigar.

With that information do the following:
  * Filter the alignment records and only take those with a mapping quality >= 30.
    (Take a look at the tutorial \ref sequence_file_section_fun_with_ranges for a reminder about using views on files)
  * For the resulting alignments, print which read was mapped against which reference id and
    the number of `seqan3::gap`s in each sequence (aligned reference and read sequence).

\hint
Given your reference sequences, you need to know which reference sequence to use for the alignment.
For this purpose, you can access `record.reference_id()` which gives you the position of the respective reference
sequence.
\endhint

Your program should print the following:

```
r001 mapped against 0 with 1 gaps in the read sequence and 2 gaps in the reference sequence.
r003 mapped against 0 with 0 gaps in the read sequence and 0 gaps in the reference sequence.
r004 mapped against 1 with 14 gaps in the read sequence and 0 gaps in the reference sequence.
```

\endassignment

\solution

\snippet doc/tutorial/10_sam_file/sam_file_solution2.cpp solution

\endsolution

# Writing alignment files

## Writing records

When writing a SAM file without any further specifications, the default file assumes that all fields are provided.
Since those are quite a lot for alignment files, we usually want to write only a subset of the data stored in the SAM format and default the rest.

For this purpose, you can use the seqan3::sam_record to write out a partial record.

\snippet doc/tutorial/10_sam_file/sam_file_writing.cpp main

Note that this only works because in the SAM format **all fields are optional**.
If we provide fewer fields when writing, default values are written.

\assignment{Assignment 3: Writing id and sequence information}

Create a small program that writes the following unmapped (see seqan3::sam_flag) read ids and sequences:

```
read1: ACGATCGACTAGCTACGATCAGCTAGCAG
read2: AGAAAGAGCGAGGCTATTTTAGCGAGTTA
```

Your ids can be of type `std::string` and your sequences of type `std::vector<seqan3::dna4>`.

Your resulting SAM file should look like this:

```
read1   4       *       0       0       *       *       0       0       ACGATCGACTAGCTACGATCAGCTAGCAG   *
read2   4       *       0       0       *       *       0       0       AGAAAGAGCGAGGCTATTTTAGCGAGTTA   *
```

\endassignment
\solution

\snippet doc/tutorial/10_sam_file/sam_file_solution3.cpp solution

\endsolution
