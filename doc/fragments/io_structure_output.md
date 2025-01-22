<!-- SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
     SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
     SPDX-License-Identifier: CC-BY-4.0
-->

### Writing structure files

Structured sequence files contain intra-molecular interactions of RNA or protein. Usually, but not necessarily, they
contain the nucleotide or amino acid sequences and descriptions as well. Interactions can be represented
either as fixed _secondary structure_, where every character is assigned at most one interaction partner
(structure of minimum free energy), or an _annotated sequence_, where every character is assigned a set
of interaction partners with specific base pair probabilities.
<br><br>
The structured sequence file abstraction supports writing ten different fields:

1. seqan3::field::seq (sequence)
2. seqan3::field::id (identifier)
3. seqan3::field::bpp (annotated sequence)
4. seqan3::field::structure (secondary structure)
5. seqan3::field::structured_seq (sequence and structure in one range)
6. seqan3::field::energy (minimum free energy)
7. seqan3::field::react (reactivity)
8. seqan3::field::react_err (reactivity error)
9. seqan3::field::comment (free text)
10. seqan3::field::offset (index of first sequence character)

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

\include test/snippet/io/structure_file/structure_file_output_temp_param_deduc.cpp

Writing to std::cout:

\include test/snippet/io/structure_file/structure_file_output_write_std_out.cpp

Note that this is not the same as writing `structure_file_output<>` (with angle brackets). In the latter case they are
explicitly set to their default values, in the former case
[automatic deduction](https://en.cppreference.com/w/cpp/language/class_template_argument_deduction) happens which
chooses different parameters depending on the constructor arguments. Prefer deduction over explicit defaults.

#### Writing record-wise

You can iterate over this file record-wise:

\include test/snippet/io/structure_file/structure_file_output_iter_by_rec.cpp

The easiest way to write to a sequence file is to use the push_back() or emplace_back() member functions. These
work similarly to how they work on a std::vector. If you pass a tuple to push_back() or give arguments to
emplace_back() the seqan3::field ID of the i-th tuple-element/argument is assumed to be the i-th value of
selected_field_ids, i.e. by default the first is assumed to be seqan3::field::seq, the second seqan3::field::id
and the third one seqan3::field::structure. You may give less fields than are selected, if the actual format you are
writing to can cope with less
(e.g. for Vienna it is sufficient to write seqan3::field::seq, seqan3::field::id and seqan3::field::structure,
even if selected_field_ids also contains seqan3::field::energy).
<br><br>
You may also use the output file's iterator for writing, however, this rarely provides an advantage.

#### Writing record-wise (custom fields)

If you want to pass a combined object for SEQ and STRUCTURE fields to push_back() / emplace_back(), or if you want
to change the order of the parameters, you can pass a non-empty fields trait object to the
structure_file_output constructor to select the fields that are used for interpreting the arguments.
<br><br>
The following snippets demonstrates the usage of such a fields trait object.

\include test/snippet/io/structure_file/structure_file_output_write_fields.cpp

A different way of passing custom fields to the file is to pass a seqan3::record – instead of a tuple – to
push_back(). The seqan3::record clearly indicates which of its elements has which seqan3::field ID so the file will
use that information instead of the template argument. This is especially handy when reading from one file and
writing to another, because you don't have to configure the output file to match the input file, it will just work:

\include test/snippet/io/structure_file/structure_file_output_pass_rec.cpp

#### Writing record-wise in batches

You can write multiple records at once, by assigning to the file:

\include test/snippet/io/structure_file/structure_file_output_mult_rec.cpp

#### File I/O pipelines

Record-wise writing in batches also works for writing from input files directly to output files, because input
files are also input ranges in SeqAn. This can be combined with file-based views to create I/O pipelines:

\include test/snippet/io/structure_file/structure_file_output_pipeline.cpp

#### Column-based writing

The record-based interface treats the file as a range of tuples (the records), but in certain situations
you might have the data as columns, i.e. a tuple-of-ranges, instead of range-of-tuples.
<br><br>
You can use column-based writing in that case, it uses operator=() and seqan3::views::zip():

\include test/snippet/io/structure_file/structure_file_output_col_based.cpp

#### Formats

Currently, the only implemented format is seqan3::format_vienna. More formats will follow soon.
