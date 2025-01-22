<!-- SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
     SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
     SPDX-License-Identifier: CC-BY-4.0
-->

### Writing SAM files

#### Construction and specialisation

The seqan3::sam_file_output class comes with two constructors, one for construction from a file name
and one for construction from an existing stream and a known format. The first
one automatically picks the format based on the extension of the file name.
The second can be used if you have a non-file stream, like std::cout or
std::ostringstream, that you want to read from and/or if you cannot use
file-extension based detection, but know that your output file has a certain
format.
<br><br>
In most cases the template parameters are deduced completely automatically:

\include test/snippet/io/sam_file/sam_file_output_filename_construction.cpp

Writing to std::cout:

\include test/snippet/io/sam_file/sam_file_output_cout_write.cpp

Note that this is not the same as writing `sam_file_output<>`
(with angle brackets). In the latter case they are explicitly set to their
default values, in the former case
[automatic deduction](https://en.cppreference.com/w/cpp/language/class_template_argument_deduction)
happens which chooses different parameters depending on the constructor arguments.
For opening from file, `sam_file_output<>` would have also worked, but for
opening from stream it would not have.

#### Writing record-wise

\include test/snippet/io/sam_file/record_based_writing.cpp

The easiest way to write to a SAM/BAM file is to use the `push_back()` member functions. These
work similarly to how they work on a std::vector.
You may also use a tuple like interface or the `emplace_back()`
function but this is not recommended since one would have to keep track of the
correct order of many fields (14 in total). For the record based interface
using `push_back()` please also see the seqan3::record documentation on how to specify
a record with the correct field and type lists.
<br><br>
You may also use the output file's iterator for writing, however, this rarely provides an advantage.

#### Writing record-wise (custom fields)

If you want to omit non-required parameter or
change the order of the parameters, you can pass a non-empty fields trait object to the
seqan3::sam_file_output constructor to select the fields that are used for interpreting the arguments.
<br><br>
The following snippet demonstrates the usage of such a `field_traits` object.

\include test/snippet/io/sam_file/record_based_writing2.cpp

A different way of passing custom fields to the file is to pass a seqan3::record – instead of a tuple – to
`push_back()`. The seqan3::record clearly indicates which of its elements has which seqan3::field so **the file will
use that information instead of the template argument**. This is especially handy when reading from one file and
writing to another, because you don't have to configure the output file to match the input file, it will just work:

\include test/snippet/io/sam_file/sam_file_output_custom_fields.cpp

This will copy the \ref seqan3::field "seqan3::field::flag" and \ref seqan3::field "seqan3::field::ref_offset" value
into the new output file.

\note Note that the other SAM columns in the output file will have a default value, so unless you specify to read
all SAM columns (see seqan3::format_sam) the output file will not be equal to the input file.

#### Writing record-wise in batches

You can write multiple records at once, by assigning to the file:

\include test/snippet/io/sam_file/sam_file_output_write_range.cpp

#### File I/O pipelines

Record-wise writing in batches also works for writing from input files directly to output files, because input
files are also input ranges in SeqAn:

\include test/snippet/io/sam_file/sam_file_output_input_range.cpp

This can be combined with file-based views to create I/O pipelines:

\snippet test/snippet/io/sam_file/sam_file_output_io_pipeline.cpp snippet

#### Formats

We currently support writing the following formats:
* seqan3::format_sam
* seqan3::format_bam
