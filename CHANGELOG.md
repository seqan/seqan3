# Changelog {#about_changelog}

[TOC]

This changelog contains a top-level entry for each release with sections on new features, API changes and notable
bug-fixes (not all bug-fixes will be listed).
See the documentation on [API stability](https://docs.seqan.de/seqan/3-master-user/about_api.html) to learn about
when API changes are allowed.

<!--
The following API changes should be documented as such:
  * a previously experimental interface now being marked as stable
  * an interface being removed
  * syntactical changes to an interface (e.g. renaming or reordering of files, functions, parameters)
  * semantical changes to an interface (e.g. a function's result is now always one larger) [DANGEROUS!]

If possible, provide tooling that performs the changes, e.g. a shell-script.
-->

# 3.0.1

## New features

#### Alphabet

* We now provide the seqan3::dna3bs alphabet that models bisulfite-treated dna4 sequence (#1191).

#### Alignment

* The function seqan3::align_pairwise can now be parallised using the seqan3::align_cfg::parallel configuration.
* The score type used in the alignment score matrix and the result type is now configurable through a template
  argument of the seqan3::align_cfg::result configuration.

#### Argument parser

* Simplified reading file extensions from formatted files with the seqan3::input_file_validator and
  seqan3::output_file_validator (#863).
* The seqan3::value_list_validator is now constructible from a range or a parameter pack.
* Enable subcommand argument parsing ([How-to](https://docs.seqan.de/seqan/3-master-user/subcommand_arg_parse.html)).
* The seqan3::argument_parser::add_option (and add_positional_option) calls do now allow enum types when using the
  seqan3::enumeration_names customisation point (#1196).

#### Build system

* `find_package(SeqAn3)` is now case-insensitive and always populates `SEQAN3_*` variables in all upper-case.

#### Core

* We now provide seqan3::lzcnt, seqan3::tzcnt, and seqan3::popcount for bit manipulation (#1141).
* Added traits for "metaprogramming" with `seqan3::type_list` and type packs.
* We now provide the SIMD functions seqan3::upcast and seqan3::upcast_signed (#1190).

#### I/O

* We now support padded alignments in the SAM/BAM format (#1173).
* We increased our input performance using a faster iterator on the stream buffer.
* Asynchronous input (background file reading) supported via seqan3::view::async_input_buffer.
* Reading field::cigar into a vector over seqan3::cigar is supported via seqan3::alignment_file_input.
* Writing field::cigar into a vector over seqan3::cigar is supported via seqan3::alignment_file_output.

### Range

* We now provide the seqan3::dynamic_bitset (#1153).
* We now provide the seqan3::views::translate_join (#1171).
* We now provide the seqan3::views::to_simd (#1190).
* We now provide the seqan3::views::kmer_hash (#946).

#### Search

* The memory footprint of FM-indices over text collections was reduced (#1363).

### Std

* We provide a `std::to_chars` overload for floating point data types in our `seqan3/std/from_chars` header (#1160).

## API changes

* **The required version of the ranges-v3 library has increased:** We now support the versions >= 0.10.0 and < 0.11.0,
  increasing the previous requirement of >= 0.5.0 and < 0.6.0.
* **Customising for third party types has changes slightly:**
  You are only affected if you added types to `seqan3::custom::`.
  Please see [About Customisation](http://docs.seqan.de/seqan/3-master-user/about_customisation.html).
* All our concept are named in the `snake_case` style now!

#### Alphabet

* The seqan3::cigar alphabet is not an seqan3::alphabet anymore but only a seqan3::semialphabet (#1285).

#### Argument parser

* The seqan3::value_list_validator is not constructible from a std::initialiser_list any more
  (e.g. `seqan3::value_list_validator{{1,2,3}}` does **not** work, use `seqan3::value_list_validator{1,2,3}` instead).
* **Changed class signature of input/output file validators:**
  Most user code will be unaffected; to fix possible compiler errors you need to add an empty template parameter list to
  the respective instances (e.g. change `input_file_validator` to `input_file_validator<>`).
* The member type that denotes which arguments a `validator` can validate has been renamed from `value_type` to
  `option_value_type`.

#### Build system

* [find_package](https://cmake.org/cmake/help/latest/command/find_package.html#version-selection) now accepts minimum
versions (e.g. `find_package(SEQAN3 3.0.1)` requires at least seqan3 with a version of `>= 3.0.1` and `< 4.0.0`).
* The variable `SEQAN3_VERSION_STRING` defined by `find_package(SEQAN3)` was renamed to `SEQAN3_VERSION`.

#### Core

* **The `type_list` header has moved:**
  If you included `<seqan3/core/type_list.hpp>` you need to change the path to `<seqan3/core/type_list/type_list.hpp>`.

#### I/O

* **Removed the field-based in- and output interface for sequence and structure files through std::get and std::tie:**
  Output can instead be achieved with seqan3::views:zip(), for input we will implement unzip() in the future.
* The `field::flag` of SAM/BAM input and output is now an **enum** instead of a simple integer (see seqan3::sam_flag).
* Uppercase seqan3::field names are deprecated. Use the lower case field names instead. You can easily find and replace
  all occurrences by the following regex: find `field::([A-Z_]+)` replace `field::\L$1`.

* **Removed the char type from the input and output files:**
  Most user code will be unaffected; however, if you have fully specified all templates of any of the input or output
  files in your code, you need to remove the template parameter to select the char type of the stream
  (e.g. change `seqan3::sequence_file_input<traits_t, fields_t, formats_t, char>` to
  `seqan3::sequence_file_input<traits_t, fields_t, formats_t>`). Before this change, setting the char type gave the
  impression that also streams over wide characters are supported which is not the case yet.

#### Range

* **The `seqan3::concatenated_sequences::data()` function has been deprecated:**
  Use `seqan3::concatenated_sequences::raw_data()` instead.

#### Search

* **Changed class signature of (bi_)fm_index:**
  All code that relies on automatic template deduction will be unaffected. In case you specified the template parameters
  of a `seqan3::fm_index` or `seqan3::bi_fm_index` you will need to add the alphabet type as first parameter and pass a
  `seqan3::text_layout` instead of a `bool` to indicate the text layout (single, collection).
  For example, `fm_index<false> index{text}` where `text` is of type `dna4_vector` needs to be changed to
  `fm_index<dna4, text_layout::single> index{text}`.

* **The `construct()` method of the (bi_)fm_index is now private:**
  Use the constructor `seqan3::fm_index::fm_index(text_t && text)` or `seqan3::bi_fm_index::bi_fm_index(text_t && text)`
  instead.

* **The `seqan3::fm_index::char_type` member was renamed to `seqan3::fm_index::alphabet_type`**
  The same applies for the `seqan3::bi_fm_index`.

* **The `seqan3::fm_index_cursor::index_char_type` member was renamed to `seqan3::fm_index_cursor::index_alphabet_type`**
  The same applies for the `seqan3::bi_fm_index_cursor`.

## Notable Bug-fixes

* All our headers are self contained now (#1085).
* The alignment algorithm with edit distance returns the correct back coordinate now (#1093)
* Inserting and deleting gaps in an empty seqan3::gap_decorator doesn't assert anymore (#1109).
* Some fixes to edge cases in BAM file writing (#1110)
* The application name of the seqan3::argument_parser is restricted to alpha-numeric characters and `_` and `-` (#1133).
* Copying and moving the `seqan3::fm_index` and `seqan3::bi_fm_index` now work properly.
* Searching in the `seqan3::fm_index` and `seqan3::bi_fm_index` constructed from a text collection containing a single
  text now return correct result.
* The view seqan3::views::take is now sized if the underlying range is sized (#1146).
* The detection of the pthread library now works correctly on linux based systems (#1200).

* Translation table for nucleotide to amino acid translation now corrected.

# 3.0.0 ("Escala")

Initial release of SeqAn3.
