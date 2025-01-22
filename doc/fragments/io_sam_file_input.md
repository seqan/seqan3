<!-- SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
     SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
     SPDX-License-Identifier: CC-BY-4.0
-->

### Reading SAM files

#### Construction and specialisation

The seqan3::sam_file_input class comes with four constructors: One for construction from a file name, one for
construction from an existing stream and a known format and both of the former with or without additional
reference information. Constructing from a file name automatically picks the format based on the extension
of the file name. Constructing from a stream can be used if you have a non-file stream, like std::cin or
std::istringstream. It also comes in handy, if you cannot use file-extension based detection, but know that
your input file has a certain format.
<br><br>
Passing reference information, e.g.

- ref_ids: The name of the references, e.g. "chr1", "chr2", ...
- ref_sequences: The reference sequence information **in the same order as the ref_ids**.

comes in handy once you want to convert the CIGAR string, read from your file, into an actual alignment.
This will be covered in the section "Transforming the CIGAR information into an actual alignment".

In most cases the template parameters are deduced automatically:

\include test/snippet/io/sam_file/sam_file_input_construction_from_filename.cpp

Reading from a std::istringstream:

\include test/snippet/io/sam_file/sam_file_input_construction_from_stream.cpp

Note that this is not the same as writing `sam_file_input<>` (with angle brackets). In the latter case they
are explicitly set to their default values, in the former case
[automatic deduction](https://en.cppreference.com/w/cpp/language/class_template_argument_deduction) happens which
chooses different parameters depending on the constructor arguments. For opening from file, `sam_file_input<>`
would have also worked, but for opening from stream it would not have.
<br><br>
You can define your own traits type to further customise the types used by and returned by this class, see
seqan3::sam_file_input_default_traits for more details. As mentioned above, specifying at least one
template parameter yourself means that you loose automatic deduction. The following is equivalent to the automatic
type deduction example with a stream from above:

\include test/snippet/io/sam_file/sam_file_input_construction_without_automatic_type_deduction.cpp

#### Reading record-wise

You can iterate over this file record-wise:

\include test/snippet/io/sam_file/sam_file_input_reading_range_based_for_loop.cpp

In the above example, `rec` has the type seqan3::sam_file_input::record_type which is a specialisation of seqan3::record
and behaves like a std::tuple (that's why we can access it via `get`). Instead of using the seqan3::field based
interface on the record, you could also use `std::get<0>` or even `std::get<dna4_vector>` to retrieve the sequence,
but it is not recommended, because it is more error-prone.

\note It is important to write `auto &` and not just `auto`, otherwise you will copy the record on every iteration.
Since the buffer gets "refilled" on every iteration, you can also move the data out of the record if you want
to store it somewhere without copying:

\include test/snippet/io/sam_file/sam_file_input_reading_move_record.cpp

#### Reading record-wise (custom fields)

If you want to skip specific fields from the record you can pass a non-empty fields trait object to the
seqan3::sam_file_input constructor to select the fields that should be read from the input. For example,
you may only be interested in the mapping flag and mapping quality of your SAM data to get some statistics.
The following snippets demonstrate the usage of such a fields trait object.

\include test/snippet/io/sam_file/sam_file_input_reading_custom_fields.cpp

When reading a file, all fields not present in the file (but requested implicitly or via the `selected_field_ids`
parameter) are ignored and the respective value in the record stays empty.

#### Reading record-wise (decomposed records)

Instead of using `get` on the record, you can also use
[structured bindings](https://en.cppreference.com/w/cpp/language/structured_binding)
to decompose the record into its elements. Considering the example of reading only the flag and mapping quality
like before you can also write:

\include test/snippet/io/sam_file/sam_file_input_reading_structured_bindings.cpp

In this case you immediately get the two elements of the tuple: `flag` of seqan3::sam_file_input::flag_type and `mapq`
of seqan3::sam_file_input::mapq_type.

\note But beware: with structured bindings you do need to get the order of elements correctly!

#### Transforming the CIGAR information into an actual alignment

In SeqAn, we represent an alignment as a tuple of two `seqan3::aligned_sequence`s.

The conversion from a CIGAR string to an alignment can be done with the function `seqan3::alignment_from_cigar`.
You need to pass the reference sequence with the position the read was aligned to and the read sequence. All
of it is already in the `record` when reading a SAM file:

\include snippet/alignment/cigar_conversion/alignment_from_cigar_io.cpp

The code will print the following:

\include snippet/alignment/cigar_conversion/alignment_from_cigar_io.err

#### Views on files

Since SeqAn files are ranges, you can also create views over files. A useful example is to filter the records
based on certain criteria, e.g. minimum length of the sequence field:

\include test/snippet/io/sam_file/sam_file_input_reading_filter.cpp

#### End of file

You can check whether a file is at its end by comparing `begin()` and `end()` (if they are the same, the file is
at its end).

#### Formats

We currently support reading the following formats:
* seqan3::format_sam
* seqan3::format_bam
