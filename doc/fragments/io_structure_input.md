<!-- SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
     SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
     SPDX-License-Identifier: CC-BY-4.0
-->

### Reading structure files

Structured sequence files contain intra-molecular interactions of RNA or protein. Usually, but not necessarily, they
contain the nucleotide or amino acid sequences and descriptions as well. Interactions can be represented
either as fixed _secondary structure_, where every character is assigned at most one interaction partner
(structure of minimum free energy), or an _annotated sequence_, where every character is assigned a set
of interaction partners with specific base pair probabilities.
<br><br>
The structured sequence file abstraction supports reading ten different fields:

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

The first three fields are retrieved by default (and in that order). The seqan3::field::structured_seq may be
selected to have sequence and structure directly stored in a more memory-efficient combined container.
If you select this field you must not select seqan3::field::seq or seqan3::field::structure.

#### Construction and specialisation

This class comes with two constructors, one for construction from a file name and one for construction from
an existing stream and a known format. The first one automatically picks the format based on the extension
of the file name. The second can be used if you have a non-file stream, like std::cin or std::istringstream,
that you want to read from and/or if you cannot use file-extension based detection, but know that your input
file has a certain format.
<br><br>
In most cases the template parameters are deduced completely automatically, e.g. reading from a std::istringstream:

\include test/snippet/io/structure_file/structure_file_input_auto_temp_deduc.cpp

Note that this is not the same as writing `structure_file_input<>` (with angle brackets). In the latter case they are
explicitly set to their default values, in the former case
[automatic deduction](https://en.cppreference.com/w/cpp/language/class_template_argument_deduction) happens which
chooses different parameters depending on the constructor arguments. For opening from file, `structure_file_input<>`
would have also worked, but for opening from stream it would not have.
<br><br>
In some cases, you do need to specify the arguments, e.g. if you want to read amino acids:

\include test/snippet/io/structure_file/structure_file_input_arg_spec.cpp

You can define your own traits type to further customise the types used by and returned by this class, see
seqan3::structure_file_default_traits_rna for more details. As mentioned above, specifying at least one
template parameter yourself means that you loose automatic deduction so if you want to read amino acids **and**
want to read from a string stream you need to give all types yourself:

\include test/snippet/io/structure_file/structure_file_input_trait_def.cpp

#### Reading record-wise

You can iterate over this file record-wise:

\include test/snippet/io/structure_file/structure_file_input_record_iter.cpp

In the above example, rec has the type seqan3::structure_file_input::record_type which is a specialisation of
seqan3::record and behaves like a std::tuple (that's why we can access it via get). Instead of using the seqan3::field
based interface on the record, you could also use `std::get<0>` or even `std::get<rna4_vector>` to retrieve the
sequence, but it is not recommended, because it is more error-prone.
<br><br>
*Note:* It is important to write `auto &` and not just `auto`, otherwise you will copy the record on every iteration.
Since the buffer gets "refilled" on every iteration, you can also move the data out of the record if you want
to store it somewhere without copying:

\include test/snippet/io/structure_file/structure_file_input_data_out.cpp

#### Reading record-wise (decomposed records)

Instead of using `get` on the record, you can also use
[structured bindings](https://en.cppreference.com/w/cpp/language/structured_binding)
to decompose the record into its elements:

\include test/snippet/io/structure_file/structure_file_input_structured_bindings.cpp

In this case you immediately get the three elements of the tuple: `seq` of seqan3::structure_file_input::seq_type,
`id` of seqan3::structure_file_input::id_type and `structure` of seqan3::structure_file_input::structure_type.
**But beware: with structured bindings you do need to get the order of elements correctly!**

#### Reading record-wise (custom fields)

If you want to skip specific fields from the record you can pass a non-empty fields trait object to the
structure_file_input constructor to select the fields that should be read from the input. For example to choose a
combined field for SEQ and STRUCTURE (see above). Or to never actually read the STRUCTURE, if you don't need it.
The following snippets demonstrate the usage of such a fields trait object.

\include test/snippet/io/structure_file/structure_file_input_skip_fields.cpp

When reading a file, all fields not present in the file (but requested implicitly or via the `selected_field_ids`
parameter) are ignored.

#### Views on files

Since SeqAn files are ranges, you can also create views over files. A useful example is to filter the records
based on certain criteria, e.g. minimum length of the sequence field:

\include test/snippet/io/structure_file/structure_file_input_filter_criteria.cpp

#### End of file

You can check whether a file is at end by comparing begin() and end() (if they are the same, the file is at end).

#### Formats

Currently, the only implemented format is seqan3::format_vienna. More formats will follow soon.
