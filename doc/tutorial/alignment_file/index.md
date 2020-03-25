# Alignment Input and Output in SeqAn {#tutorial_alignment_file}

***Learning Objective:***

<b>Learning Objective:</b> <br/>
You will get an overview of how to read and write alignment files.
This tutorial is a walk-through with links into the API documentation and also meant as a source for copy-and-paste code.

\tutorial_head{Medium, 60 min, \ref setup\, \ref tutorial_alphabets\, \ref tutorial_sequence_file,}

[TOC]

# Introduction

Alignment files are used to store pairwise alignments between two (biological) sequences.
Common file formats are the Sequence Alignment/Map format (SAM) and BLAST output format.
Next to the alignment, those formats store additional information like the start positions or mapping qualities.
Alignment files are a little more complex than sequence files but the basic design is the same.
If you are new to SeqAn, we strongly recommend to do the tutorial \ref tutorial_sequence_file first.

# Alignment file formats

## SAM format

SAM stands for Sequence Alignment/Map format. It is a TAB-delimited text format consisting of a header
section, which is optional, and an alignment section
(see the [official SAM specifications](https://samtools.github.io/hts-specs/SAMv1.pdf)).

Here is an example of a SAM file:

\include doc/tutorial/alignment_file/example.sam

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

BAM is the binary format version of SAM. It provides the same data as the SAM format
with the only difference that **the header is mandatory**.

# Alignment file fields

\copydetails seqan3::alignment_file_input::field_ids

Note that some of the fields are specific to the SAM format, while some are specific to BLAST.
To make things clearer, here is the table of SAM columns connected to the corresponding alignment file field:

| #  | SAM Column ID |  FIELD name                                                                       |
|:--:|:--------------|:----------------------------------------------------------------------------------|
| 1  | QNAME         | seqan3::field::id                                                                 |
| 2  | FLAG          | seqan3::field::flag                                                               |
| 3  | RNAME         | seqan3::field::ref_id                                                             |
| 4  | POS           | seqan3::field::ref_offset                                                         |
| 5  | MAPQ          | seqan3::field::mapq                                                               |
| 6  | CIGAR         | implicitly stored in seqan3::field::alignment or directly in seqan3::field::cigar |
| 7  | RNEXT         | seqan3::field::mate (tuple pos 0)                                                 |
| 8  | PNEXT         | seqan3::field::mate (tuple pos 1)                                                 |
| 9  | TLEN          | seqan3::field::mate (tuple pos 2)                                                 |
| 10 | SEQ           | seqan3::field::seq                                                                |
| 11 | QUAL          | seqan3::field::qual                                                               |

## File extensions

The formerly introduced formats can be identified by the following file name extensions
(this is important for automatic format detection from a file name as you will learn in the next section).

| File Format  | File Extensions   |
| -------------|-------------------|
| SAM          |   `.sam`          |
| BAM          |   `.bam`          |

You can access and modify the valid file extensions via the `file_extension` member variable in a format tag:

\snippet doc/tutorial/alignment_file/alignment_file_snippets.cpp file_extensions

# Reading alignment files

Before we start, you should copy and paste this [example file](example.sam) into a file location of your choice
(we use `/tmp/` in the examples, so make sure you adjust your path).

\attention Make sure the file you copied is tab delimited!

## Construction

The construction works analogously to sequence files by passing a file name,
in which case all template parameters are automatically deduced (by the file name extension).
Or you can pass a stream (e.g. std::cin or std::stringstream), but then you need to know your format beforehand:

\snippet doc/tutorial/alignment_file/alignment_file_snippets.cpp main
\snippet doc/tutorial/alignment_file/alignment_file_snippets.cpp filename_construction
\snippet doc/tutorial/alignment_file/alignment_file_snippets.cpp main_end

## Reading custom fields

In many cases you are not interested in all of the information in a file.
For this purpose, we provide the possibility to select specific seqan3::field's for a file.
The file will read only those fields and fill the record accordingly.

You can select fields by providing a seqan3::fields object as an extra parameter to the constructor:

\snippet doc/tutorial/alignment_file/alignment_file_snippets.cpp main
\snippet doc/tutorial/alignment_file/alignment_file_snippets.cpp read_custom_fields
\snippet doc/tutorial/alignment_file/alignment_file_snippets.cpp main_end

\attention The order in which you specify the selected fields determines the order of elements in the seqan3::record.

In the example above we only select the id, sequence and flag information
so the seqan3::record object has three tuple elements that are decomposed using structural bindings.

Note that this is possible for all SeqAn file objects.

\assignment{Assignment 1: Accumulating mapping qualities}

Let's assume we want to compute the average mapping quality of a SAM file.

For this purpose, write a small program that
    * only reads the mapping quality (field::mapq) out of a SAM file and
    * computes the average of all qualities.

Use the following file to test your program:

\snippet doc/tutorial/alignment_file/alignment_file_solution1.cpp sam_file

It should output:
```
Average: 27.4
```

\endassignment
\solution

\snippet doc/tutorial/alignment_file/alignment_file_solution1.cpp solution

\endsolution

# Alignment representation in the SAM format

In SeqAn we represent an alignment as a tuple of two *aligned_sequences*,
as you have probably learned by now from the alignment tutorial.

The SAM format is the common output format of read mappers where you align short read sequences
to one or more large reference sequences.
In fact, the SAM format stores those alignment information only partially:
It **does not store the reference sequence** but only the read sequence and a *CIGAR* string
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
Fields set to `0` or `*` are defaulted and contain no information.

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

By default, the `seqan3::alignment_file_input` will always read the `seqan3::field::cigar` and store it
into a `std::vector<seqan3::cigar>`:

\snippet doc/tutorial/alignment_file/alignment_file_read_cigar.cpp code

## Reading the CIGAR information into an actual alignment

In SeqAn, the conversion from a CIGAR string to an alignment (two *aligned_sequences*) is done automatically for you. You can access it by querying `seqan3::field::alignment` from the record:

\snippet doc/tutorial/alignment_file/alignment_file_snippets.cpp main
\snippet doc/tutorial/alignment_file/alignment_file_snippets.cpp alignments_without_ref
\snippet doc/tutorial/alignment_file/alignment_file_snippets.cpp main_end

In the example above, you can only safely access the aligned read.

\attention The **unknown** aligned reference sequence at the first position in the alignment tuple cannot be accessed
           (e.g. via the `operator[]`). It is represented by a dummy type that throws on access.

Although the SAM format does not handle reference sequence information,
you can provide these information to the seqan3::alignment_file_input which automatically fills the alignment object.
You can pass reference ids and reference sequences as additional constructor parameters:

\snippet doc/tutorial/alignment_file/alignment_file_snippets.cpp main
\snippet doc/tutorial/alignment_file/alignment_file_snippets.cpp alignments_with_ref
\snippet doc/tutorial/alignment_file/alignment_file_snippets.cpp main_end

\assignment{Assignment 2: Combining sequence and alignment files}

Read in the following reference sequence FASTA file (see the sequence file tutorial if you need a remainder):

\snippet doc/tutorial/alignment_file/alignment_file_solution2.cpp ref_file

Then read in the following SAM file while providing the reference sequence information.
Only read in the id, reference id, mapping quality, and alignment.

\snippet doc/tutorial/alignment_file/alignment_file_solution2.cpp sam_file

With those information do the following:
  * Filter the alignment records and only take those with a mapping quality >= 30.
    (Take a look at the tutorial \ref sequence_file_section_fun_with_ranges for a reminder how to use views on files)
  * For the resulting alignments, print which read was mapped against with reference id and
    the number of seqan3::gap's in each sequence (aligned reference and read sequence).

\note reference ids (field::ref_id) are given as an index of type `std::optional<int32_t>`
      that denote the position of the reference id in the `ref_ids` vector passed to the alignment file.

Your program should print the following:

```
r001 mapped against 0 with 1 gaps in the read sequence and 2 gaps in the reference sequence.
r003 mapped against 0 with 0 gaps in the read sequence and 0 gaps in the reference sequence.
r004 mapped against 1 with 14 gaps in the read sequence and 0 gaps in the reference sequence.
```

\endassignment

\solution

\snippet doc/tutorial/alignment_file/alignment_file_solution2.cpp solution

\endsolution

# Writing alignment files

## Writing records

When writing a SAM file without any further specifications, the default file assumes that all fields are provided.
Since those are quite a lot for alignment files, we usually want to write only a subset of the data stored in the SAM format and default the rest.

For this purpose, you can also select specific fields by giving an additional seqan3::fields object to the constructor.

\attention The **order** of the field tags in your seqan3::fields object will determine
           the order of values stored in the record type!

\snippet doc/tutorial/alignment_file/alignment_file_snippets.cpp main
\snippet doc/tutorial/alignment_file/alignment_file_snippets.cpp writing
\snippet doc/tutorial/alignment_file/alignment_file_snippets.cpp main_end

Note that this only works because in the SAM format **all fields are optional**.
So if we provide less fields when writing, default values are printed.

\assignment{Assignment 3: Writing id and sequence information}

Write a small program that writes the following read ids + sequences:

```
read1: ACGATCGACTAGCTACGATCAGCTAGCAG
read2: AGAAAGAGCGAGGCTATTTTAGCGAGTTA

```

Your ids can be of type `std::string` and your sequences of type `std::vector<seqan3::dna4>`.

Your resulting SAM file should look like this:

```
read1   0       *       0       0       *       *       0       0       ACGATCGACTAGCTACGATCAGCTAGCAG   *
read2   0       *       0       0       *       *       0       0       AGAAAGAGCGAGGCTATTTTAGCGAGTTA   *
```

\endassignment
\solution

\include doc/tutorial/alignment_file/alignment_file_solution3.cpp

\endsolution
