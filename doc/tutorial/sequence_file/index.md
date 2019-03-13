# Sequence File Input and Output {#tutorial_sequence_file}

<b>Learning Objective:</b> <br/>
You will get an overview of how file Input/Output is handled in seqan3 and learn how to read and write
sequence files. This tutorial is a walk-through with links into the API documentation and also meant as a
source for copy-and-paste code.

\tutorial_head{Easy, 30 min, \ref setup\, alphabets\, ranges, [POSIX conventions](https://www.math.uni-hamburg.de/doc/java/tutorial/essential/attributes/_posix.html)}

[TOC]

# File I/O in SeqAn3

Most file formats in bioinformatics are structured as lists of records.
In SeqAn3 we model our files as a **range over records**.
This interface allows us to easily stream over a file, apply filters and convert formats,
sometimes only in a single line of code.
The file format is automatically detected by the file name extension and compressed files can be handled
without any effort. We can even stream over files in a python-like way with range-based for loops:

```cpp
for (auto & record : file)
    // do something with my record
```

We will explain the details about reading and writing files in the [Sequence File](#section_sequence_files) section below.
Currently, SeqAn3 supports the following file formats:

- Sequence file formats:
  - seqan3::sequence_file_format_fasta
  - seqan3::sequence_file_format_fastq

- Structure file formats:
  - seqan3::structure_file_format_vienna

- Alignment file formats:
  - seqan3::alignment_file_format_sam

\warning Access to compressed files relies on external libraries.
For instance, you need to have *zlib* installed for reading `.gz` files and *libbz2* for reading `.bz2` files.
You can check whether you have installed these libraries by running `cmake .` in your build directory.
If `-- Optional dependency: ZLIB-x.x.x found.` is displayed on the command line then you can read/write
compressed files in your programs. TODO what about bz2

## Basic layout of SeqAn3 file objects

Before we dive into the details, we will outline the general design of our file objects,
hoping that it will make the following tutorial easier to understand.

As mentioned above, our file object is a range over records.
More specifically over objects of type seqan3::record which is basically just a std::tuple that holds the data.
To identify or specialise which data is read/written and contained in the records,
we use seqan3::field tags (e.g. seqan3::field::SEQ denotes sequence information).
The seqan3::field tags are shared between file formats and allow for easy file conversion.

Output files can handle various types that fulfill the requirements of the format (e.g.
a sequence has to be a range over an alphabet).
In contrast to this, input files have certain default types for record fields that can be modified via
a *traits type*. For example, on construction you can specify seqan3::sequence_file_default_traits_dna
or seqan3::sequence_file_default_traits_aa for reading `dna` and `protein` sequences respectively (section
[File traits](#section_file_traits) will covers this in more detail).

Opening and closing files is also handled automatically.
If a file cannot be opened for reading or writing, a seqan3::file_open_error is thrown.

# Sequence Files {#section_sequence_files}

Sequence files are the most generic and common biological files.
Well-known formats include FASTA and FASTQ,
but some may also be interested in treating SAM or BAM files as sequence files, discarding the alignment.
The Sequence file abstraction supports reading four different fields:

  1. seqan3::field::SEQ
  2. seqan3::field::ID
  3. seqan3::field::QUAL
  4. seqan3::field::SEQ_QUAL (sequence and qualities in one range, see seqan3::qualified for details)

The first three fields are retrieved by default (and in that order!).
The last field may be selected to have sequence and qualities directly stored in a more memory-efficient
combined container.
If you select the last field you may not select seqan3::field::SEQ or seqan3::field::QUAL.

# Sequence file formats

### FASTA format

A FASTA record contains the sequence id and the sequence characters. Here is an example of a FASTA file:

```
>seq1
CCCCCCCCCCCCCCC
>seq2
CGATCGATC
```
### FASTQ format

A FASTQ record contains an additional quality value for each sequence character. Here is an example of a FASTQ file:

```
@seq1
CCCCCCCCCCCCCCC
+
IIIIIHIIIIIIIII
@seq2
CGATCGATC
+
IIIIIIIII
```

### File extensions

The formerly introduced formats can be identified by the following file name extensions
(this is important for automatic format detection from a file name as you will learn in the next section).

| File Format  | File Extensions   |
| -------------|-------------------|
| FASTA        |   `.fa`, `.fasta`, `.fna`, `.ffn`, `.ffa`, `.frn` |
| FASTQ        |   `.fq`, `.fastq` |

You can access the valid file extension via the `file_extension` member variable in a format:

\snippet doc/tutorial/sequence_file/sequence_file_snippets.cpp file_extensions

You can also customise this list if you want to allow different or additional file extensions:

\snippet doc/tutorial/sequence_file/sequence_file_snippets.cpp modify_file_extensions

# Reading a sequence file

You can include the SeqAn3 sequence file functionality with:

\snippet doc/tutorial/sequence_file/sequence_file_snippets.cpp include

## Construction

At first, you need to construct a seqan3::sequence_file_input object that handles file access.
In most cases you construct from a file name:

\snippet test/snippet/io/sequence_file/sequence_file_input.cpp template_deduction

All template parameters of the seqan3::sequence_file_input are automatically deduced, even the format!
We **detect the format by the file name extension**.
The file extension in the example above is `.fasta` so we choose the seqan3::sequence_file_format_fasta.

You can also construct a sequence file object directly from a stream (e.g. std::cin or std::stringstream),
but then you need to know your format beforehand:

\snippet doc/tutorial/sequence_file/sequence_file_snippets.cpp construct_from_cin

### File traits {#section_file_traits}

The seqan3::sequence_file_input needs to know the types of the data you are reading in
(e.g. that your sequence is a `dna` sequence) at **compile-time**.
These necessary types are defined in the  **traits type** argument (seqan3::sequence_file_input::traits_type)
which by default is set to seqan3::sequence_file_input_default_traits_dna.
We thereby assume that you want to read the sequence into a std::vector over a seqan3::dna5 alphabet
and store the sequence names/id in a std::string.

In case you want to read a **protein** sequence instead we also provide the
seqan3::sequence_file_input_default_traits_aa traits type which sets the SEQ field to a std::vector over
the seqan3::aa27 alphabet.

You can specify the traits object as the first template argument for the sequence file:

\snippet doc/tutorial/sequence_file/sequence_file_snippets.cpp amino_acid_type_trait

You can also customise the types by inheriting from one of the default traits and changing the type manually.
See the detailed information on seqan3::sequence_file_input_default_traits_dna for an example.

## Reading records

After construction you can now read the sequence records. As described in the basic file layout,
our file objects behave like ranges so you can use a range based for loop to conveniently iterate over the file:

\snippet test/snippet/io/sequence_file/sequence_file_input.cpp record_iter

In the above example, `rec` has the type seqan3::sequence_file_input::record_type
which is a specialisation of seqan3::record and behaves like an std::tuple
(that's why we can access it via `get`).

\note It is important to write `auto &` and not just `auto`, otherwise you will copy the record on every iteration.
      Since the buffer gets "refilled" on every iteration, you can also move the data out of the record if you want
      to store it somewhere without copying:
      \snippet test/snippet/io/sequence_file/sequence_file_input.cpp auto_ref


Since the return type seqan3::record behaves like a tuple, you can also use
[structured bindings](http://en.cppreference.com/w/cpp/language/structured_binding)
to decompose the record into its elements:

\snippet test/snippet/io/sequence_file/sequence_file_input.cpp decomposed

In this case you immediately get the two elements of the tuple:
`seq` of seqan3::sequence_file_input::sequence_type and `id` of seqan3::sequence_file_input::id_type.
**But beware: with structured bindings you do need to get the order of elements correctly!**

You can read up more on the different ways to stream over the file object in the detailed documentation
of seqan3::sequence_file_input.

\assignment{Exercise: Reading a FASTQ file}
Copy and paste the following FASTQ file to some location, e.g. the tmp directory:

\snippet doc/tutorial/sequence_file/sequence_file_snippets.cpp fastq_file

and then create a simple program that reads in all the records from that file and prints the id,
sequence and quality information to the command line using the seqan3::debug_stream.
Do not use the structured bindings for now but access the record via `get<>()`!

Note: you include the seqan3::debug_stream with the following header

\snippet doc/tutorial/sequence_file/sequence_file_snippets.cpp include_debug_stream

\endassignment
\solution

\snippet doc/tutorial/sequence_file/sequence_file_solution1.cpp solution

The code will print the following:
```bash
ID:  seq1
SEQ: AGCTAGCAGCGATCG
QUAL: IIIIIHIIIIIIIII
ID:  seq2
SEQ: CGATCGATC
QUAL: IIIIIIIII
ID:  seq3
SEQ: AGCGATCGAGGAATATAT
QUAL: IIIIHHGIIIIHHGIIIH
```
\endsolution

## Reading only a subset of information

In some cases, you might be only interested in a subset of your data,
e.g. you only need the sequence information but not the id and quality.
For this purpose, you can specialise the seqan3::sequence_file_input::selected_fields type by giving an
additional seqan3::fields object to the constructor.

\attention The **order** of the field tags in your seqan3::fields object will determine
           the order of values stored in the record type!

\snippet doc/tutorial/sequence_file/sequence_file_snippets.cpp custom_fields

<a name="assignment_reading_seq_qual"></a>
\assignment{Exercise: Reading a FASTQ file into a SEQ_QUAL object}

While working on the same FASTQ file from assignment 1,
create a new program that reads in the `ID` and the `SEQ_QUAL` object and
prints the sequence and quality information to the command line.
Take a look at the seqan3::qualified type on how to handle the combination of sequence and quality information.

\hint
Pipe the seqan3::view::get on the `SEQ_QUAL` object to access the sequence or quality information .
\endhint

\endassignment
\solution

\snippet doc/tutorial/sequence_file/sequence_file_solution2.cpp solution

The code will print the following:
```bash
ID:  seq1
SEQ: AGCTAGCAGCGATCG
QUAL: IIIIIHIIIIIIIII
ID:  seq2
SEQ: CGATCGATC
QUAL: IIIIIIIII
ID:  seq3
SEQ: AGCGATCGAGGAATATAT
QUAL: IIIIHHGIIIIHHGIIIH
```
\endsolution

# Sequence files as views

Since SeqAn files are ranges, you can also create views over files.
This enables us to create solutions for lot of use cases using only a few lines of code.

## Reading a file in batches

A common use case is to read chunks from a file instead of the whole file at once.

You can do so easily on a file range by using the seqan3::view::take.

\snippet doc/tutorial/sequence_file/sequence_file_snippets.cpp read_in_batches

The example above will iterate over the file, by reading 2 records at a time.
If no two records are available any more, it will just print the remaining single record.

## Applying a filter to a file

In some occasions you are only interested in sequence records that fulfill a certain criteria,
e.g. having a minimum sequence length or a minimum quality.
Just like in the example with the *take* view you can use std::ranges::filter for this purpose:

\snippet test/snippet/io/sequence_file/sequence_file_input.cpp file_view

To remind you what you have learned in the ranges & views tutorial before,
a view is not applied immediately but **lazy evaluated**.
That means your file is still parsed record by record and not at once.

\assignment{Exercise: Fun with ranges}

Implement a small program that reads in a FASTQ file and outputs the first 2 sequences
that have an average quality of at least 40.

Hints:
* You may want to (but do not need to) use `std::accumulate` when aggregating values.
* You need the following includes for view::filter
\snippet doc/tutorial/sequence_file/sequence_file_snippets.cpp include_filter
* You need the following includes for ranges::view::take
\snippet doc/tutorial/sequence_file/sequence_file_snippets.cpp include_ranges_take

Test your program on the following FASTQ file:

\snippet doc/tutorial/sequence_file/sequence_file_solution3.cpp fastq_file

It should output:
```console
seq1
seq3
```
\endassignment
\solution
\snippet doc/tutorial/sequence_file/sequence_file_solution3.cpp solution
\endsolution

# Writing a sequence file

## Construction

You construct the seqan3::sequence_file_output just like the seqan3::sequence_file_input
by giving it a file name or a stream.

\snippet test/snippet/io/sequence_file/sequence_file_output.cpp template_deduction

Writing to std::cout:

\snippet test/snippet/io/sequence_file/sequence_file_output.cpp cout_write

## Writing records

The easiest way to write to a sequence file is to use the seqan3::sequence_file_output::push_back()
or seqan3::sequence_file_output::emplace_back() member functions.
These work similarly to how they work on an std::vector.

\snippet test/snippet/io/sequence_file/sequence_file_output.cpp record_wise_iteration

If you pass a tuple to `push_back()` or give arguments to `emplace_back()` the order of elements is assumed
to be the same as the one in the seqan3::sequence_file_output::selected_field_ids.
For the above example the default FASTA fields are first seqan3::field::SEQ,
second seqan3::field::ID and the third one seqan3::field::QUAL.
You may give less fields than are selected if the actual format you are writing to can cope with less
(e.g. for FastA it is sufficient to give sequence and name information).

Of course you can customise the seqan3::sequence_file_output::selected_field_ids as well.

\snippet doc/tutorial/sequence_file/sequence_file_snippets.cpp writing_custom_fields

\assignment{Exercise: Naive quality trimming}

Implement a naive quality trimming approach
that reads in a FASTQ file and writes out another FASTQ file with the trimmed reads.
The quality trimming should scan the reads from left to right
and cut off the sequence whenever the quality drops below 10 (phred score).

Test your code on the following file:

\snippet doc/tutorial/sequence_file/sequence_file_solution4.cpp fastq_file

It should create a new file that looks like this:

```
@seq1
CGATCG
+
IIIIII
@seq2
AGCGATCGAG
+
IIIIHHGIII
@seq3
AGCTAGCA
+
IIIIIHII
```

\endassignment
\solution
\snippet doc/tutorial/sequence_file/sequence_file_solution4.cpp solution
\endsolution

## Files as views

Again we want to point out a convenient advantages of modelling files as ranges.
In the "reading a file" section you already saw a few examples of how to pipe a view into
a seqan3::sequence_file_input object. In the same way you can pipe the output file:

\snippet doc/tutorial/sequence_file/sequence_file_snippets.cpp piping_in_out

\assignment{Exercise: Naive quality trimming 2 }

Copy and paste the solution of the previous exercise and remove the for loop in favour of a pipe notation.

The result should be the same.

\endassignment
\solution
\snippet doc/tutorial/sequence_file/sequence_file_solution5.cpp solution
\endsolution

# File conversion

As mentioned before, the seqan3::field tags are shared between formats which allows for easy file conversion.
For example you can read in a FASTQ file and output a FASTA file in one line:

\snippet doc/tutorial/sequence_file/sequence_file_snippets.cpp file_conversion

Yes that's it! Of course this only works because all fields that are required in FASTA are provided in FASTQ.
The other way around would not work as easily because we have no quality information (and would make less sense too).
