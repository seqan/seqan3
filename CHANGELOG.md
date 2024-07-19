# Changelog {#about_changelog}

<!--
    SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
    SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
    SPDX-License-Identifier: CC-BY-4.0
-->

[TOC]

This changelog contains a top-level entry for each release with sections on new features, API changes and notable
bug-fixes (not all bug-fixes will be listed).

Get to know SeqAn3 with our [tutorials](https://docs.seqan.de/seqan3/main_user/usergroup1.html).

Please see the release announcement: https://www.seqan.de/announcing-seqan3/

See the porting guide for some help on porting: https://docs.seqan.de/seqan3/main_user/howto_porting.html

See the documentation on [API stability](https://docs.seqan.de/seqan3/main_user/about_api.html) to learn about
when API changes are allowed.

<!--
The following API changes should be documented as such:
  * a previously experimental interface now being marked as stable
  * an interface being removed
  * syntactical changes to an interface (e.g. renaming or reordering of files, functions, parameters)
  * semantic changes to an interface (e.g. a function's result is now always one larger) [DANGEROUS!]

If possible, provide tooling that performs the changes, e.g. a shell-script.
-->

# 3.3.1

## New features

## Notable Bug-fixes

## API changes

#### Dependencies
  * We now use Doxygen version 1.9.8 to build our documentation ([\#3197](https://github.com/seqan/seqan3/pull/3197)).
  * We bumped the minimal CMake version to 3.5 ([\#3223](https://github.com/seqan/seqan3/pull/3223)).

# 3.3.0

## New features

#### Alignment
  * The function `seqan3::alignment_from_cigar` creates an alignment (tuple of two aligned sequences) from a
    CIGAR vector (`std::vector<seqan3::cigar>`) ([\#3057](https://github.com/seqan/seqan3/pull/3057)) or a
    CIGAR string (`std::string`) ([\#3077](https://github.com/seqan/seqan3/pull/3077)).
  * The function `seqan3::cigar_from_alignment` creates a CIGAR vector (`std::vector<seqan3::cigar>`) from an alignment
    (tuple of 2 aligned sequences) ([\#3057](https://github.com/seqan/seqan3/pull/3057)).

#### Alphabet
  * Improved performance of vector assignment for alphabets ([\#3038](https://github.com/seqan/seqan3/pull/3038)).
  * Improved performance of `seqan3::dna4::complement()` ([\#3026](https://github.com/seqan/seqan3/pull/3026)).
  * Char literals returning `std::vector` are now `constexpr` if supported by the compiler
    ([\#3073](https://github.com/seqan/seqan3/pull/3073)).

#### I/O
  * Made `seqan3::sam_file_header::program_info_t` easier to copy ([\#3145](https://github.com/seqan/seqan3/pull/3145)).

#### Search

  * Added a constructor to the `seqan3::interleaved_bloom_filter` for decompressing a compressed
    `seqan3::interleaved_bloom_filter` ([\#3082](https://github.com/seqan/seqan3/pull/3082)).

## Notable Bug-fixes

#### Alignment
  * Resolved an issue resulting in a wrong alignment score ([\#3098](https://github.com/seqan/seqan3/pull/3098)).

#### I/O
  * Empty SAM/BAM files must at least write a header to ensure a valid file
    ([\#3081](https://github.com/seqan/seqan3/pull/3081)).
  * Reading SAM/BAM files is 2x faster than before
    ([\#3106](https://github.com/seqan/seqan3/pull/3106)).
  * `seqan3::sequence_file_input_traits` now allows `char` as a sequence alphabet
    ([\#3128](https://github.com/seqan/seqan3/pull/3128)).

#### Utility
  * Fixed spin delay having no effect on the PowerPC platform ([\#3129](https://github.com/seqan/seqan3/pull/3129)).

## API changes

#### Alignment
  * The fields `seqan3::field::offset` and `seqan3::field::alignment` have been removed from `seqan3::sam_record`
    ([\#3058](https://github.com/seqan/seqan3/pull/3058), [\#3089](https://github.com/seqan/seqan3/pull/3089)).
    For `seqan3::field::offset`, please check the soft clipping of the CIGAR string (`seqan3::sam_record::cigar()`).
    For `seqan3::field::alignment`, please use `seqan3::alignment_from_cigar`.

#### Compiler
  * Dropped support for GCC 10 ([\#3148](https://github.com/seqan/seqan3/pull/3148)).
  * Added support for GCC 13 ([\#3148](https://github.com/seqan/seqan3/pull/3148)).

#### Dependencies
  * We require at least CMake 3.16 for our test suite. Note that the minimum requirement for using SeqAn3 is unchanged
    ([\#3050](https://github.com/seqan/seqan3/pull/3050)).
  * We now use Doxygen version 1.9.6 to build our documentation ([\#3116](https://github.com/seqan/seqan3/pull/3116)).
  * Updated sdsl-lite to 3.0.3 ([\#3170](https://github.com/seqan/seqan3/pull/3170),
    [\#3174](https://github.com/seqan/seqan3/pull/3174)).
  * Compatibility with SeqAn2: SeqAn2's namespace was changed from `seqan` to `seqan2`. For interoperability, an
    up-to-date checkout of [SeqAn2's main branch](https://github.com/seqan/seqan) is required
    ([\#3156](https://github.com/seqan/seqan3/pull/3156)).

# 3.2.0

## New features

#### Alphabet
  * `seqan3::cigar` can now be assigned from `std::string_view` ([\#2966](https://github.com/seqan/seqan3/pull/2966)).
  * Added `seqan3::views::char_strictly_to`. Behaves like `seqan3::views::char_to`, but throws on invalid
    input ([\#2898](https://github.com/seqan/seqan3/pull/2898)).

#### I/O
 * Added `seqan3::sequence_file_option::fasta_ignore_blanks_before_id` to ignore blanks before IDs when reading FASTA
   files. E.g., `>       some_id` will only store `"some_id"` as ID
   ([\#2770](https://github.com/seqan/seqan3/pull/2770)).

#### Search
  * Improved performance of `seqan3::counting_vector::operator+=` ([\#2930](https://github.com/seqan/seqan3/pull/2930)).

#### Utility
  * Added `seqan3::list_traits::repeat` ([\#2899](https://github.com/seqan/seqan3/pull/2899)).

## Notable Bug-fixes

#### Core
  * Added missing implementations for AVX512 ([\#2920](https://github.com/seqan/seqan3/pull/2920) and
    [\#2926](https://github.com/seqan/seqan3/pull/2926)).

#### IO
  * FASTA files containing IDs starting with `>`, e.g., `> >MyID`, are now parsed correctly
    ([\#2869](https://github.com/seqan/seqan3/pull/2869)).

#### Search
  * Relaxed `kmer_hash_view::iterator` difference requirement ([\#2931](https://github.com/seqan/seqan3/pull/2931)).
  * Relaxed `seqan3::views::minimiser` requirements to be C++20-compatible
    ([\#2845](https://github.com/seqan/seqan3/pull/2845)).
  * Relaxed `seqan3::views::kmer_hash` requirements to be C++20-compatible
    ([\#2843](https://github.com/seqan/seqan3/pull/2843)).

#### Utility

  * `seqan3::views::single_pass_input` cannot propagate the `std::ranges::output_range` property because it cannot
    satisfy the following requirement ([\#2775](https://github.com/seqan/seqan3/pull/2775)):
    ```cpp
    *it++ = value;
    // must be the same as
    *it = value; ++it;
    // but it actually would be the same as
    ++it; *it = value;
    ```
  * Fixed signature of `seqan3::detail::store_sse4`. This might have affected some public API
    ([\#2893](https://github.com/seqan/seqan3/pull/2893)).
  * Relaxed `seqan3::views::to_simd` requirements to be C++20-compatible
    ([\#2849](https://github.com/seqan/seqan3/pull/2849)).

## API changes
GCC 7, 8, and 9 have been removed. All headers in `seqan3/std/` except `charconv` and `new` have been deprecated, please
use the equivalent `std` includes.
The namespace `std::cpp20` has been deprecated, please use `std::`.

`seqan3::views::to` has been changed to `seqan3::ranges::to`. Since it is not a view anymore, it cannot be properly
deprecated. Please keep this in mind if you encounter errors with `seqan3::views::to`.

#### Compiler
  * Dropped support for GCC 7 and 8 ([\#2891](https://github.com/seqan/seqan3/pull/2891)).
  * Dropped support for GCC 9 ([\#2952](https://github.com/seqan/seqan3/pull/2952)).
  * Removed C++17 support ([\#2915](https://github.com/seqan/seqan3/pull/2915)).

#### I/O
 * Changed the default of `output_options::fasta_blank_before_id` to `false`
   ([\#2769](https://github.com/seqan/seqan3/pull/2769)).
 * Changed the default of `bgzf_thread_count` to `4`. Previously, all available threads on the machine were utilised
   ([\#2911](https://github.com/seqan/seqan3/pull/2911)).
 * The global variable `bgzf_thread_count` is not static anymore. Modifying the variable now affects the runtime of the
   entire program. Formerly, this variable was specific to each translation unit
   ([\#2752](https://github.com/seqan/seqan3/pull/2752)).

#### Range
  * Replaced `seqan3::views::to` (implemented via range-v3) with `seqan3::ranges::to` (implemented in SeqAn3).
    `seqan3::ranges::to` provides a subset of C++23's `std::ranges::to` and will be replaced with the STL-equivalent
    in a future version ([\#2969](https://github.com/seqan/seqan3/pull/2969)).
  * Replaced the implementation of `seqan3::views::chunk`. It is now implemented in SeqAn3 and does not use
    the range-v3 implementation anymore. `seqan3::views::chunk` is similar to C++23's `std::views::chunk`
    and will be replaced with the STL-equivalent in a future version
    ([\#2975](https://github.com/seqan/seqan3/pull/2975)).
  * Replaced the implementation of `seqan3::views::join_with`. It is now implemented in SeqAn3 and does not use
    the range-v3 implementation anymore. `seqan3::views::join_with` is equivalent to C++23's `std::views::join_with`
    and will be replaced with the STL-equivalent in a future version
    ([\#2973](https://github.com/seqan/seqan3/pull/2973)).
  * Replaced the implementation of `seqan3::views::zip`. It is now implemented in SeqAn3 and does not use
    the range-v3 implementation anymore. `seqan3::views::zip` is equivalent to C++23's `std::views::zip`
    and will be replaced with the STL-equivalent in a future version
    ([\#2971](https://github.com/seqan/seqan3/pull/2971)).

#### Dependencies
  * We now use Doxygen version 1.9.4 to build our documentation ([\#2979](https://github.com/seqan/seqan3/pull/2979)).
  * Removed range-v3 ([\#2998](https://github.com/seqan/seqan3/pull/2998)).
  * Updated cereal to 1.3.2 ([\#3012](https://github.com/seqan/seqan3/pull/3012)).
  * Updated sdsl-lite to 3.0.1 ([\#3012](https://github.com/seqan/seqan3/pull/3012)).

# 3.1.0

## New features

#### Build system

* We now use Doxygen version 1.9.2 to build our documentation ([\#2765](https://github.com/seqan/seqan3/pull/2765)).

## Notable Bug-fixes

#### Argument Parser

* Fixed an issue that led to the wrong option type being printed on errors
  ([\#2836](https://github.com/seqan/seqan3/pull/2836)).

#### Search

* Resolved an issue that prevented the FM-Index from being constructed correctly for inputs larger than 4 GiB
  ([\#2756](https://github.com/seqan/seqan3/pull/2756)).

## API changes

The files deprecated in [3.0.3](#api303) (denoted by `[deleted without replacement]`) have been removed.

# 3.0.3

Note that 3.1.0 will be the first API stable release and interfaces in this release might still change.

* Check out our updated [SeqAn3 Cookbook](https://docs.seqan.de/seqan/3.0.3/cookbook.html). It contains a listing of code
  examples on how to perform particular tasks using the library.

* SeqAn 3.0.3 is known to compile with GCC 7.5, 8.4, 9.3, 10.3, and 11.1. Future versions (e.g. GCC 11.2 and 12) might
  work, but were not yet available at the time of this release.

## New features

#### Alphabet

* Added `seqan3::phred94`, a quality type that represents the full Phred Score range (Sanger format) and is used for
  PacBio Phred scores of HiFi reads ([\#2290](https://github.com/seqan/seqan3/pull/2290)).

#### Argument Parser

* The `seqan3::argument_parser` constructor was changed to accept a strong-type `seqan3::update_notifications::(on|off)`
  instead of a simple bool (which was subject to unexpected implicit conversion)
  ([\#2180](https://github.com/seqan/seqan3/pull/2180)).
* We expanded the `seqan3::output_file_validator`, with a parameter `seqan3::output_file_open_options` to allow overwriting
  output files ([\#2009](https://github.com/seqan/seqan3/pull/2009)).
* The `seqan3::argument_parser` has a new member function `seqan3::argument_parser::is_option_set` that
  checks whether an option, identified by its long or short name, was set on the command line by the user
  ([\#1859](https://github.com/seqan/seqan3/pull/1859)).
* The subcommand of the `seqan3::argument_parser` may contain a dash
  ([\#2306](https://github.com/seqan/seqan3/pull/2306)).

#### Build system

* We now use Doxygen version 1.9.1 to build our documentation ([\#2327](https://github.com/seqan/seqan3/pull/2327)).

#### I/O

* Explicit record-classes with explicit member accessor for our file implementations. We added `seqan3::sequence_record`
  for `seqan3::sequence_file_(in|out)put`, `seqan3::sam_record` for `seqan3::sam_file_(in|out)put` and
  `seqan3::structure_record` for `seqan3::structure_file_(in|out)put`. You can now access the `id` in a sequence file
  (e.g. `fasta` file) record via `record.id()` instead of `seqan3::get<seqan3::field::id>(record)`. This will allow us
  to add convenient functions that compute information based on the record itself and to provide better documentation.
  ([\#2340](https://github.com/seqan/seqan3/pull/2340), [\#2380](https://github.com/seqan/seqan3/pull/2380),
  [\#2389](https://github.com/seqan/seqan3/pull/2389))

#### Search

* The `seqan3::fm_index_cursor` and `seqan3::bi_fm_index_cursor` can be serialised
  ([\#2048](https://github.com/seqan/seqan3/pull/2048)).
* The `seqan3::fm_index_cursor` exposes its suffix array interval ([\#2076](https://github.com/seqan/seqan3/pull/2076)).
* The `seqan3::interleaved_bloom_filter` supports counting occurrences of a range of values
  ([\#2373](https://github.com/seqan/seqan3/pull/2373)).
* The `seqan3::interleaved_bloom_filter` supports clearing of bins
  ([\#2428](https://github.com/seqan/seqan3/pull/2428)).

## Notable Bug-fixes

#### Argument Parser

* Argument parsing with subcommands: If the user forgets or misspells the subcommand, the error is thrown when calling
  `seqan3::argument_parser::parse()` and not on construction of the `seqan3::argument_parser`
  ([\#2179](https://github.com/seqan/seqan3/pull/2179)).
* The `seqan3::regex_validator` parses `std::filesystem::path`'s correctly now
  ([\#2216](https://github.com/seqan/seqan3/pull/2216)).
* The `seqan3::input_file_validator` and `seqan3::input_file_validator` support extensions containing a dot
  ([\#2363](https://github.com/seqan/seqan3/pull/2363)).
* The Argument Parser accepts containers of all values it is able to parse, e.g. a `std::vector` of enums or `bool`
  ([\#2381](https://github.com/seqan/seqan3/pull/2381)).
* The Argument Parser's help pages now contain author and email information
  ([\#2374](https://github.com/seqan/seqan3/pull/2374)).
* The Argument Parser is able to process defaults for list arguments
  ([\#2394](https://github.com/seqan/seqan3/pull/2394)).

#### Range

* The `seqan3::views::kmer_hash` does not return wrong values when combined with `std::views::reverse` on a text of the
  same size as the kmer ([\#2416](https://github.com/seqan/seqan3/pull/2416)).

#### I/O

* The SAM/BAM format reads the quality field (or any other text field) that starts with an asterisk (`*`) but is not
  empty correctly now ([\#2184](https://github.com/seqan/seqan3/pull/2184)).
* Requesting the alignment without also requesting the sequence for BAM files containing empty CIGAR strings does now
  not result in erroneous parsing ([\#2418](https://github.com/seqan/seqan3/pull/2418)).
* BAM files with 64 references are now parsed correctly ([\#2423](https://github.com/seqan/seqan3/pull/2423)).
* BAM files not containing a plain text header are now accepted ([\#2536](https://github.com/seqan/seqan3/pull/2536)).
* Writing `gz`-compressed output no longer results in `bgzf`-compressed output. This change may have following effects
  ([\#2458](https://github.com/seqan/seqan3/pull/2458)):
  * A noticeable slowdown when writing `gz`-compressed content since, in contrast to `bgzf`, `gz` does not feature
    parallelisation (magnitude depends on the application and level of parallelisation).
  * A reduced output size when writing `gz`-compressed content due to `gz` storing less metadata than `bgzf`
    (up to 20% smaller file size).
  * The processed data should experience no negative effects since `gz` and `bgzf` are **fully compatible**.
  * This bug may also cause unexpected parallelisation when reading `gz`-compressed input. This is the case when the
    `gz`-compressed input was also generated with SeqAn.
* The following requirements of parsing a SAM/BAM header where relaxed as they were in violation of the specification
  ([\#2605](https://github.com/seqan/seqan3/pull/2605)):
  * The order of tags within the header may be arbitrary. Before, necessary tags where required to be first.
  * Unknown tags will emit a warning. Before, a error was raised.

## API changes {#api303}

Most of our API or header file changes will trigger a deprecation warning to let you know if something changed and, if
applicable, when it will be removed. We recommend upgrading version-by-version to check whether you need to change code.
You can either directly check the reported code or verify with our documentation how the new API should be used.

For a complete list of behavioural changes in our public and internal API, you can consult our API stability
regression test suite and patches at https://github.com/seqan/seqan3/tree/main/test/api_stability/3.0.2.

#### Alignment

* The member constants of `seqan3::aminoacid_similarity_matrix` were changed to lower case
  ([\#2599](https://github.com/seqan/seqan3/pull/2599)):
  * `seqan3::aminoacid_similarity_matrix::BLOSUM30` is replaced by `seqan3::aminoacid_similarity_matrix::blosum30`.
  * `seqan3::aminoacid_similarity_matrix::BLOSUM45` is replaced by `seqan3::aminoacid_similarity_matrix::blosum45`.
  * `seqan3::aminoacid_similarity_matrix::BLOSUM62` is replaced by `seqan3::aminoacid_similarity_matrix::blosum62`.
  * `seqan3::aminoacid_similarity_matrix::BLOSUM80` is replaced by `seqan3::aminoacid_similarity_matrix::blosum80`.
* Deprecated seqan3::alignment_coordinate and it will be removed in 3.1.0
  ([\#2602](https://github.com/seqan/seqan3/pull/2602)).

Header Changes:

```cpp
#include <seqan3/alignment/band/static_band.hpp> [deleted without replacement]
#include <seqan3/alignment/matrix/advanceable_alignment_coordinate.hpp> [deleted without replacement]
#include <seqan3/alignment/scoring/gap_scheme.hpp> [deleted without replacement]
```
#### Alphabet

* We relaxed some requirements of `seqan3::alphabet_base<alphabet_t>`
  ([\#2427](https://github.com/seqan/seqan3/pull/2427)):
  * Old requirements: `alphabet_t::rank_to_char` and `alphabet_t::char_to_rank` must be lookup tables.
  * New requirements: `alphabet_t::rank_to_char` and `alphabet_t::char_to_rank` must be static member functions.

  This allows for more flexible rank <-> char conversion implementations. Lookup tables are still possible within those
  static member functions. However, alphabets that do not need a lookup table can now use easier and/or more efficient
  implementations. For example, `seqan3::gap` always returns rank `0` or char `-`, or `seqan3::phred42` where the rank
  and char representations are offset by a fixed value.
* We relaxed a requirement of `seqan3::nucleotide_base<alphabet_t>`
  ([\#2584](https://github.com/seqan/seqan3/pull/2584)):
  * Old requirement: `alphabet_t::complement_table` must be a lookup table.
  * New requirement: `alphabet_t::rank_complement` must be a static member function.
* Removed seqan3::char_is_valid_for requirement from seqan3::writable_alphabet and
  seqan3::detail::writable_constexpr_alphabet
  ([\#2337](https://github.com/seqan/seqan3/pull/2337)).
* Removed seqan3::cigar_op, please use seqan3::cigar::operation instead
  ([\#2388](https://github.com/seqan/seqan3/pull/2388)).
* The literal 'M'_cigar_op was renamed to 'M'_cigar_operation
  ([\#2388](https://github.com/seqan/seqan3/pull/2388)).
* Renamed `seqan3::phred68legacy` to `seqan3::phred68solexa` ([\#2522](https://github.com/seqan/seqan3/pull/2522)).
* Renamed `seqan3::sam_dna16` to `seqan3::dna16sam`
  ([\#2521](https://github.com/seqan/seqan3/pull/2521)).
* Removed phred constructors ([\#2537](https://github.com/seqan/seqan3/pull/2537)):
  * Use `seqan3::phred42::assign_phred()` or `seqan3::operator""_phred42` instead of `seqan3::phred42(phred_type)`.
  * Use `seqan3::phred63::assign_phred()` or `seqan3::operator""_phred63` instead of `seqan3::phred63(phred_type)`.
  * Use `seqan3::phred94::assign_phred()` or `seqan3::operator""_phred94` instead of `seqan3::phred94(phred_type)`.
  * Use `seqan3::phred68legacy::assign_phred()` or `seqan3::operator""_phred68legacy` instead of
    `seqan3::phred68legacy(phred_type)`.
* Renamed `seqan3::quality_base` to `seqan3::phred_base`
  ([\#2539](https://github.com/seqan/seqan3/pull/2539)).
* Added the `seqan3::literals` namespace containing all literals. This adds the option to use
  `using namespace seqan3::literals` to import literal operators. The old way of explicitly importing specific
  operators via `using seqan3::operator""_{dna4, rna4, ...}` is not affected by this change
  ([\#2568](https://github.com/seqan/seqan3/pull/2568)).
 * **Breaking change**: Swapped the meaning of `seqan3::alphabet_variant::is_alternative` and
   `seqan3::alphabet_variant::holds_alternative` ([\#2596](https://github.com/seqan/seqan3/pull/2596)).

Header Changes:

```cpp
#include <seqan3/alphabet/cigar/{cigar_op => cigar}.hpp>
#include <seqan3/alphabet/nucleotide/{sam_dna16 => dna16sam}.hpp}>
#include <seqan3/alphabet/quality/{phred68legacy => phred68solexa}.hpp>
#include <seqan3/alphabet/quality/{quality_base => phred_base}.hpp>
```

#### Argument Parser

* `seqan3::output_file_validator` cannot be constructed with the extension list alone anymore, you need to specify one
  of the seqan3::output_file_open_options options.
  ([\#2009](https://github.com/seqan/seqan3/pull/2009)).
* The enum names of `seqan3::option_spec` were changed to lower case
  ([\#2285](https://github.com/seqan/seqan3/pull/2285)):
  * `seqan3::option_spec::DEFAULT` is replaced by `seqan3::option_spec::standard`.
  * `seqan3::option_spec::REQUIRED` is replaced by `seqan3::option_spec::required`.
  * `seqan3::option_spec::ADVANCED` is replaced by `seqan3::option_spec::advanced`.
  * `seqan3::option_spec::HIDDEN` is replaced by `seqan3::option_spec::hidden`.

#### Core

* Deprecated seqan3::range_compatible_concept and it will be removed in 3.1.0
  ([\#2265](https://github.com/seqan/seqan3/pull/2265)).

Header Changes:

```cpp
#include <seqan3/core/algorithm/bound.hpp> [Functionality included in alignment/configuration/align_config_band.hpp]
#include <seqan3/core/{algorithm => configuration}/configuration.hpp>
#include <seqan3/core/{algorithm => configuration}/pipeable_config_element.hpp>

#include <seqan3/{core => utility}/char_operations/predicate.hpp>
#include <seqan3/{core => utility}/char_operations/transform.hpp>

#include <seqan3/{core => utility/tuple}/common_tuple.hpp>

#include <seqan3/{core/concept/tuple => utility/tuple/concept}.hpp>

#include <seqan3/{core => utility}/math.hpp>

#include <seqan3/{core => utility/tuple}/pod_tuple.hpp>

#include <seqan3/{core => utility}/simd/concept.hpp>
#include <seqan3/{core => utility}/simd/{simd_algorithm => algorithm}.hpp>
#include <seqan3/{core => utility}/simd/simd.hpp>
#include <seqan3/{core => utility}/simd/simd_traits.hpp>
#include <seqan3/{core => utility}/simd/{view_iota_simd => views/iota_simd}.hpp>
#include <seqan3/{core => utility}/simd/{view_to_simd => views/to_simd}.hpp>

#include <seqan3/{core/tuple_utility => utility/tuple/pop_front}.hpp>
#include <seqan3/{core/tuple_utility => utility/tuple/split}.hpp>

#include <seqan3/{core => utility}/type_list/traits.hpp>
#include <seqan3/{core => utility}/type_list/type_list.hpp>

#include <seqan3/{core => utility}/type_traits/basic.hpp>
#include <seqan3/{core => utility}/type_traits/concept.hpp>
#include <seqan3/{core => utility}/type_traits/function{ =>_traits}.hpp>
#include <seqan3/{core => utility}/type_traits/{lazy => lazy_conditional.hpp>
#include <seqan3/{core/type_traits/pack => utility/type_pack/traits}.hpp>
#include <seqan3/{core => utility}/type_traits/pre.hpp>  [deleted without replacement]
#include <seqan3/core/{type_traits/range => range/type_traits}.hpp>
```

#### I/O

* Deprecated `seqan3::field::seq_qual`. Use `seqan3::field::seq` and `seqan3::field::qual` instead.
  ([\#2379](https://github.com/seqan/seqan3/pull/2379)). Check out
  [SeqAn3 Cookbook - Write Record](https://docs.seqan.de/seqan/3.0.3/cookbook.html#autotoc_md51) for usage.
* Renamed seqan3::alignment_file\* to seqan3::sam_file\*
  ([\#2459](https://github.com/seqan/seqan3/pull/2459)):
  * `seqan3::alignment_file_header` is replaced by `seqan3::sam_file_header`.
  * `seqan3::alignment_file_input_default_traits` is replaced by `seqan3::sam_file_input_default_traits`.
  * `seqan3::alignment_file_input` is replaced by `seqan3::sam_file_input`.
  * `seqan3::alignment_file_input_format` is replaced by `seqan3::sam_file_input_format`.
  * `seqan3::alignment_file_input_options` is replaced by `seqan3::sam_file_input_options`.
  * `seqan3::alignment_file_output` is replaced by `seqan3::sam_file_output`.
  * `seqan3::alignment_file_output_format` is replaced by `seqan3::sam_file_output_format`.
  * `seqan3::alignment_file_output_options` is replaced by `seqan3::sam_file_output_options`.
* `seqan3::sam_file_input` and `seqan3::sam_file_output` do not accept `seqan3::field::ref_seq`,
  `seqan3::field::evalue` and `seqan3::field::bit_score` anymore.
  ([\#2658](https://github.com/seqan/seqan3/pull/2658)).
* The `seqan3::get` accessor for I/O records, e.g. `seqan3::get<seqan3::field::id>(record)`, is deprecated, please use
  the corresponding member accessor ([\#2420](https://github.com/seqan/seqan3/pull/2420)):
  * If you used files as views with `seqan3::views::get<seqan3::field::id>` to project a single field, e.g.
    * `seqan3::views::get<seqan3::field::id>(fin)` => `std::views::transform(fin, [](auto && record){ return record.id(); })`
    * `fin | seqan3::views::get<seqan3::field::id>()` => `fin | std::views::transform([](auto && record){ return record.id(); })`
    * or per projection: `fin | std::views::transform(&decltype(fin)::%record_type::id)`
  * `seqan3::sequence_record`:
    * `seqan3::get<seqan3::field::id>(record)` => `record.id()`
    * `seqan3::get<seqan3::field::seq>(record)` => `record.sequence()`
    * `seqan3::get<seqan3::field::qual>(record)` => `record.base_qualities()`
  * `seqan3::structure_record`:
    * `seqan3::get<seqan3::field::id>(record)` => `record.id()`
    * `seqan3::get<seqan3::field::seq>(record)` => `record.sequence()`
    * `seqan3::get<seqan3::field::structure>(record)` => `record.sequence_structure()`
    * `seqan3::get<seqan3::field::energy>(record)` => `record.energy()`
    * `seqan3::get<seqan3::field::bpp>(record)` => `record.base_pair_probability_matrix()`
  * `seqan3::sam_record`:
    * `seqan3::get<seqan3::field::id>(record)` => `record.id()`
    * `seqan3::get<seqan3::field::seq>(record)` => `record.sequence()`
    * `seqan3::get<seqan3::field::qual>(record)` => `record.base_qualities()`
    * `seqan3::get<seqan3::field::offset>(record)` => `record.sequence_position()`
    * `seqan3::get<seqan3::field::alignment>(record)` => `record.alignment()`
    * `seqan3::get<seqan3::field::ref_id>(record)` => `record.reference_id()`
    * `seqan3::get<seqan3::field::ref_offset>(record)` => `record.reference_position()`
    * `seqan3::get<seqan3::field::header_ptr>(record)` => `record.header_ptr()`
    * `seqan3::get<seqan3::field::flag>(record)` => `record.flag()`
    * `std::get<0>(seqan3::get<seqan3::field::mate>(record))` => `record.mate_reference_id()`
    * `std::get<1>(seqan3::get<seqan3::field::mate>(record))` => `record.mate_position()`
    * `std::get<2>(seqan3::get<seqan3::field::mate>(record))` => `record.template_length()`
    * `seqan3::get<seqan3::field::mapq>(record)` => `record.mapping_quality()`
    * `seqan3::get<seqan3::field::cigar>(record)` => `record.cigar_sequence()`
    * `seqan3::get<seqan3::field::tags>(record)` => `record.tags()`

Header Changes:

```cpp
#include <seqan3/io/{alignment_file=> sam_file}/format_bam.hpp>
#include <seqan3/io/{alignment_file=> sam_file}/format_sam.hpp>
#include <seqan3/io/{alignment_file=> sam_file}/header.hpp>
#include <seqan3/io/{alignment_file=> sam_file}/input.hpp>
#include <seqan3/io/{alignment_file=> sam_file}/input_format_concept.hpp>
#include <seqan3/io/{alignment_file=> sam_file}/input_options.hpp>
#include <seqan3/io/{alignment_file/misc => sam_file/sam_flag}.hpp>
#include <seqan3/io/{alignment_file=> sam_file}/output.hpp>
#include <seqan3/io/{alignment_file=> sam_file}/output_format_concept.hpp>
#include <seqan3/io/{alignment_file=> sam_file}/output_options.hpp>
#include <seqan3/io/{alignment_file=> sam_file}/sam_tag_dictionary.hpp>
```
#### Range

* We made `seqan3::views::convert` NOAPI and moved it to `seqan3/utility/views/convert.hpp`. You can still use
  `seqan3::views::convert` in the meantime, but we encourage using `std::views::transform` instead as shown in our
  [Cookbook](https://docs.seqan.de/seqan/3.0.3/cookbook.html#cookbook_convert_alphabet_range)
  ([\#2524](https://github.com/seqan/seqan3/pull/2524)).
* Deprecated `seqan3::views::as_const`, there is no alternative other than reimplementing it yourself
  ([\#2567](https://github.com/seqan/seqan3/pull/2567)).
* Deprecated `seqan3::views::drop`, use `std::views::drop` or `seqan3::views::type_reduce | std::views::drop`.
  ([\#2540](https://github.com/seqan/seqan3/pull/2540))
* Deprecated `seqan3::views::join`. Please use `std::views::join` or `seqan3::views::join_with` instead
  ([\#2526](https://github.com/seqan/seqan3/pull/2526)).
* Deprecated `seqan3::views::move`, use the `std::ranges::move` algorithm, `std::[cpp20::]move_iterator` or an explicit
  for loop where you move the value.
  ([\#2563](https://github.com/seqan/seqan3/pull/2563))
* Deprecated `seqan3::views::take` and it will be removed in 3.1.0. Use `std::views::take` instead
  ([\#2541](https://github.com/seqan/seqan3/pull/2541)).
* Deprecated `seqan3::views::take_line` and it will be removed in 3.1.0
  ([\#2525](https://github.com/seqan/seqan3/pull/2525)).
* Deprecated `seqan3::views::take_exactly`. Please use `std::views::take` or `std::views::counted` instead
  ([\#2601](https://github.com/seqan/seqan3/pull/2601)).
* Deprecated `seqan3::views::take_until` and it will be removed in 3.1.0. Use
  `std::views::take_while(std::not_fn(predicate))` instead ([\#2604](https://github.com/seqan/seqan3/pull/2604)).
* Deprecated `seqan3::views::take_until_and_consume` and it will be removed in 3.1.0. There is no alternative other
  than reimplementing it yourself ([\#2604](https://github.com/seqan/seqan3/pull/2604)).
* Deprecated `seqan3::views::to_upper` and it will be removed in 3.1.0, use
  `std::views::transform([](auto && chr){return std::toupper(chr)})`.
  ([\#2538](https://github.com/seqan/seqan3/pull/2538))
* Deprecated `seqan3::views::to_lower` and it will be removed in 3.1.0, use
   `std::views::transform([](auto && chr){return std::tolower(chr)})`.
   ([\#2556](https://github.com/seqan/seqan3/pull/2556))
* Deprecated `seqan3::views::persist`. There is no replacement, use lvalues instead of rvalues
  ([\#2553](https://github.com/seqan/seqan3/pull/2553)).
* Renamed `seqan3::gap_decorator::unaligned_seq_type` to `seqan3::gap_decorator::unaligned_sequence_type`
  ([\#2564](https://github.com/seqan/seqan3/pull/2564)).
* Renamed `seqan3::views::get` to `seqan3::views::elements` ([\#2554](https://github.com/seqan/seqan3/pull/2554)).
* Renamed `seqan3::translation_frames::*FRAME*` ([\#2565](https://github.com/seqan/seqan3/pull/2565)):
  * `seqan3::translation_frames::FWD_FRAME_0` is replaced by `seqan3::translation_frames::forward_frame0`.
  * `seqan3::translation_frames::FWD_FRAME_1` is replaced by `seqan3::translation_frames::forward_frame1`.
  * `seqan3::translation_frames::FWD_FRAME_2` is replaced by `seqan3::translation_frames::forward_frame2`.
  * `seqan3::translation_frames::REV_FRAME_0` is replaced by `seqan3::translation_frames::reverse_frame0`.
  * `seqan3::translation_frames::REV_FRAME_1` is replaced by `seqan3::translation_frames::reverse_frame1`.
  * `seqan3::translation_frames::REV_FRAME_2` is replaced by `seqan3::translation_frames::reverse_frame2`.
  * `seqan3::translation_frames::FWD_REV_0` is replaced by `seqan3::translation_frames::forward_reverse0`.
  * `seqan3::translation_frames::FWD_REV_1` is replaced by `seqan3::translation_frames::forward_reverse1`.
  * `seqan3::translation_frames::FWD_REV_2` is replaced by `seqan3::translation_frames::forward_reverse2`.
  * `seqan3::translation_frames::FWD` is replaced by `seqan3::translation_frames::forward_frames`.
  * `seqan3::translation_frames::REV` is replaced by `seqan3::translation_frames::reverse_frames`.
  * `seqan3::translation_frames::SIX_FRAME` is replaced by `seqan3::translation_frames::six_frames`.
* Renamed `seqan3::type_reduce_view` to `seqan3::type_reduce_t` ([\#2587](https://github.com/seqan/seqan3/pull/2587)).

Header Changes:

```cpp
#include <seqan3/{range/concept => alphabet/range/sequence}.hpp>
#include <seqan3/{range/concept => alphabet/range/concept}.hpp>

#include <seqan3/{range => utility}/container/aligned_allocator.hpp>
#include <seqan3/{range => alphabet}/container/{bitcompressed_vector => bitpacked_sequence}.hpp>
#include <seqan3/{range => alphabet}/container/concatenated_sequences.hpp>
#include <seqan3/{range => utility}/container/concept.hpp>
#include <seqan3/{range => utility}/container/dynamic_bitset.hpp>
#include <seqan3/{range => utility}/container/small_string.hpp>
#include <seqan3/{range => utility}/container/small_vector.hpp>

#include <seqan3/{range => alignment}/decorator/gap_decorator.hpp>

#include <seqan3/{range => alphabet/range}/hash.hpp>

#include <seqan3/{range => io}/views/async_input_buffer.hpp>
#include <seqan3/{range => alphabet}/views/char_to.hpp>
#include <seqan3/{range => utility}/views/chunk.hpp>
#include <seqan3/{range => alphabet}/views/complement.hpp>
#include <seqan3/{range => utility}/views/convert.hpp>
#include <seqan3/{range => utility}/views/deep.hpp>
#include <seqan3/{range => utility}/views/enforce_random_access.hpp>
#include <seqan3/{range => utility}/views/{get => elements}.hpp>
#include <seqan3/{range => utility}/views/interleave.hpp>
#include <seqan3/range/views/istreambuf.hpp> [deleted without replacement]
#include <seqan3/{range => utility}/views/{join => join_with}.hpp>
#include <seqan3/{range => search}/views/kmer_hash.hpp>
#include <seqan3/{range => search}/views/minimiser.hpp>
#include <seqan3/{range => search}/views/minimiser_hash.hpp>
#include <seqan3/{range => utility}/views/pairwise_combine.hpp>
#include <seqan3/range/views/persist> [deleted without replacement]
#include <seqan3/{range => alphabet}/views/rank_to.hpp>
#include <seqan3/{range => utility}/views/repeat.hpp>
#include <seqan3/{range => utility}/views/repeat_n.hpp>
#include <seqan3/{range => utility}/views/single_pass_input.hpp>
#include <seqan3/{range => utility}/views/slice.hpp>
#include <seqan3/range/views/take.hpp> [deleted without replacement]
#include <seqan3/range/views/take_exactly.hpp> [deleted without replacement]
#include <seqan3/range/views/take_line.hpp> [deleted without replacement]
#include <seqan3/range/views/take_until.hpp> [deleted without replacement]
#include <seqan3/{range => utility}/views/to.hpp>
#include <seqan3/{range => alphabet}/views/to_char.hpp>
#include <seqan3/{range => alphabet}/views/to_rank.hpp>
#include <seqan3/{range => alphabet}/views/translate.hpp>
#include <seqan3/{range => alphabet}/views/translate_join.hpp>
#include <seqan3/{range => alphabet}/views/trim_quality.hpp>
#include <seqan3/{range => utility}/views/type_reduce.hpp>
#include <seqan3/{range => utility}/views/zip.hpp>
```

#### Search

* We removed the concepts `seqan3::[bi_]fm_index[_cursor]_specialisation`. We did this because we currently have only one
  implementation modelling each concept and are not completely sure if the current definition of the concepts is the
  right one. If you used those concepts, you can check whether the cursor type is `seqan3::[bi_]fm_index_cursor` as a
  substitute. ([\#2348](https://github.com/seqan/seqan3/pull/2348))

# 3.0.2

Note that 3.1.0 will be the first API stable release and interfaces in this release might still change.

* Check out our new [SeqAn3 Cookbook](https://docs.seqan.de/seqan/3.0.2/cookbook.html). It contains a listing of code
  examples on how to perform particular tasks using the library.

* SeqAn 3.0.2 is known to compile with GCC 7.5, 8.4, 9.3 and 10.2. Future versions (e.g. GCC 10.3 and 11) might work,
  but were not yet available at the time of this release.

## New features

#### Alignment

* The alignment algorithm can now be invoked with a user defined callback function using the alignment configuration
  `seqan3::align_cfg::on_result` ([\#1876](https://github.com/seqan/seqan3/pull/1876)).
* The function `seqan3::align_pairwise` accepts a `std::pair` of sequences as input
  ([\#1913](https://github.com/seqan/seqan3/pull/1913)).
* We lowered the requirements of the `seqan3::aligned_sequence` concept by removing everything that needs write
  access to the object. We then added a new `seqan3::writable_aligned_sequence` concept which extends
  `seqan3::aligned_sequence` with the requirements that need write access (e.g. `insert_gap`)
  ([\#1933](https://github.com/seqan/seqan3/pull/1933)).

#### Argument Parser

* The following functions accept a `seqan3::argument_parser::option_spec::ADVANCED` to control what is
  displayed on the (advanced) help page:
  * `seqan3::argument_parser::add_section`
  * `seqan3::argument_parser::add_subsection`
  * `seqan3::argument_parser::add_line`
  * `seqan3::argument_parser::add_list_item`
  Note that other `seqan3::argument_parser::option_spec`s like `REQUIRED` are ignored
  ([\#1652](https://github.com/seqan/seqan3/pull/1652)).

#### I/O

* The `seqan3::format_fasta` accepts the file extension `.fas` as a valid extension for the FASTA format
  ([\#1599](https://github.com/seqan/seqan3/pull/1599)).

#### Build system

* Add top-level `CMakeLists.txt` ([\#1475](https://github.com/seqan/seqan3/pull/1475)).
* We now use Doxygen version 1.8.20 to build our documentation ([\#2081](https://github.com/seqan/seqan3/pull/2081)).

#### Range

* The `seqan3::views::minimiser` has been added. This is a view that computes the minimum in a window shifted over a
  range of comparable values ([\#1654](https://github.com/seqan/seqan3/pull/1654)).
* The `seqan3::views::minimiser_hash` has been added. This is a view that computes the minimisers of a range of type
  `seqan3::semialphabet` ([\#1721](https://github.com/seqan/seqan3/pull/1721)).

#### Search

* Added `seqan3::interleaved_bloom_filter`, a data structure that efficiently answers set-membership queries for
  multiple bins ([\#920](https://github.com/seqan/seqan3/pull/920)).
* Added `seqan3::search_cfg::hit`, which allows dynamic configuration of the hit strategy.
  ([\#1853](https://github.com/seqan/seqan3/pull/1853)).
* Added `seqan3::search_cfg::on_result`, which allows providing a custom callback for the search algorithm
  ([\#2019](https://github.com/seqan/seqan3/pull/2019)).

## API changes

* The required version of the ranges-v3 library has increased:
  We now support the versions >= 0.11.0 and < 0.12.0, increasing the previous requirement of >= 0.10.0 and < 0.11.0
  ([\#2014](https://github.com/seqan/seqan3/pull/2014)).

#### Alignment

* The alignment configuration elements have been refactored:
  * All options are now classes and must be constructed explicitly even if they don’t take any arguments,
    e.g. `seqan3::align_cfg::vectorise` -> `seqan3::align_cfg::vectorised{}`.
  * The configuration `seqan3::align_cfg::band` has been replaced by `seqan3::align_cfg::band_fixed_size`,
    and will directly be initialised with a `seqan3::align_cfg::lower_diagonal` and `seqan3::align_cfg::upper_diagonal`
    instead of a `seqan3::static_band` class. It also directly exposes the `lower_diagonal` and `upper_diagonal` as
    public members ([\#1873](https://github.com/seqan/seqan3/pull/1873)).
  * The configuration `seqan3::align_cfg::mode` has been replaced by two separate configuration elements
    `seqan3::align_cfg::method_global` and `seqan3::align_cfg::method_local`
    ([\#1918](https://github.com/seqan/seqan3/pull/1918)).
  * The configuration `seqan3::align_cfg::aligned_ends` has been replaced by `seqan3::align_cfg::method_global`. The
    free end-gaps are now initialised via constructor arguments to the `seqan3::align_cfg::method_global` configuration
    ([\#2119](https://github.com/seqan/seqan3/pull/2119)).
  * The configuration `seqan3::align_cfg::vectorise` has been replaced by `seqan3::align_cfg::vectorised`
    ([\#2026](https://github.com/seqan/seqan3/pull/2026)).
  * The configuration `seqan3::align_cfg::scoring` has been replaced by `seqan3::align_cfg::scoring_scheme`
    ([\#2027](https://github.com/seqan/seqan3/pull/2027)).
  * The configuration `seqan3::align_cfg::result` has been replaced by
    [`seqan3::align_cfg::output_*` options](https://docs.seqan.de/seqan/3.0.2/group__alignment.html).
    When no output configuration was configured, the default behaviour changed from computing only the score to
    all possible outputs. Please read the linked documentation above carefully to understand all implied changes
    ([\#2024](https://github.com/seqan/seqan3/pull/2024) & [\#2035](https://github.com/seqan/seqan3/pull/2035)).
  * The configuration `seqan3::align_cfg::gap` has been replaced by `seqan3::align_cfg::gap_cost_affine`, which is
    directly initialised with the relevant gap scores ([\#2037](https://github.com/seqan/seqan3/pull/2037)).
  * The configuration `seqan3::align_cfg::max_error` has been replaced by `seqan3::align_cfg::min_score`, and thus
    prepares it for non-edit scoring schemes in the future as well
    ([\#2021](https://github.com/seqan/seqan3/pull/2021)).

Header Changes:

```cpp
#include <seqan3/alignment/configuration/{align_config_max_error.hpp => align_config_min_score.hpp}>
#include <seqan3/alignment/configuration/{align_config_scoring.hpp => align_config_scoring_scheme.hpp}>
#include <seqan3/alignment/configuration/{align_config_vectorise.hpp => align_config_vectorised.hpp}>
#include <seqan3/alignment/configuration/{align_config_gap.hpp => align_config_gap_cost_affine.hpp}>
#include <seqan3/alignment/configuration/{align_config_mode.hpp => align_config_method.hpp}>
#include <seqan3/alignment/configuration/{align_config_result.hpp => align_config_score_type.hpp}>

#include <seqan3/alignment/configuration/align_config_aligned_ends.hpp> [deleted without replacement; is now part of seqan3::align_cfg::method_global]

#include <seqan3/{alignment/pairwise => core/algorithm}/alignment_range.hpp>
```

#### Core

* In accordance with the standard, the following concepts have been renamed:
  * `std::default_constructible` to `std::default_initializable`
  * `std::readable` to `std::indirectly_readable`
  * `std::writable` to `std::indirectly_writable` ([\#1860](https://github.com/seqan/seqan3/pull/1860)).
* The `seqan3::remove_cvref_t` has been replaced by `std::remove_cvref_t`
  ([\#2079](https://github.com/seqan/seqan3/pull/2079)).

#### Range

* The `seqan3::begin()`, `seqan3::end()`, `seqan3::cbegin()`, `seqan3::cend()`, `seqan3::size()`, `seqan3::empty()`
  functions have been deprecated. Use `std::ranges::{begin|end|cbegin|cend|size|empty}()` instead
  ([\#1663](https://github.com/seqan/seqan3/pull/1663)).
* The `seqan3::forward_range` has been removed. Use `std::ranges::borrowed_range` instead
  ([\#2038](https://github.com/seqan/seqan3/pull/2038)).
* The `seqan3::views:trim` has been renamed to `seqan3::views:trim_quality`
  ([\#2025](https://github.com/seqan/seqan3/pull/2025)).

Header Changes:

```cpp
#include <seqan3/range/views/{trim.hpp => trim_quality.hpp}>
```

#### Search

* Moved `seqan3::search` from `search/algorithm/` to `search/` ([\#1696](https://github.com/seqan/seqan3/pull/1696)).
* The `seqan3::search_result_range` returns now a `seqan3::search_result` which unifies the interface for all the search
  instances, e.g. using an index over a single text or a text collection
  ([\#1706](https://github.com/seqan/seqan3/pull/1706)).
* Configuration refactoring:
  * The configuration `seqan3::search_cfg::max_error` has been replaced by individual configuration elements:
    * `seqan3::search_cfg::max_error{seqan3::search_cfg::total}` to `seqan3::search_cfg::max_error_total{}`
    * `seqan3::search_cfg::max_error{seqan3::search_cfg::insertion}` to `seqan3::search_cfg::max_error_insertion{}`
    * `seqan3::search_cfg::max_error{seqan3::search_cfg::deletion}` to `seqan3::search_cfg::max_error_deletion{}`
    * `seqan3::search_cfg::max_error{seqan3::search_cfg::substitution}` to
      `seqan3::search_cfg::max_error_substitution{}` ([\#1861](https://github.com/seqan/seqan3/pull/1861)).
  * The max error configurations can be initialised with either a `seqan3::search_cfg::error_rate` or
    `seqan3::search_cfg::error_count`, and can be reassigned ([\#1861](https://github.com/seqan/seqan3/pull/1861)).
  * The configuration `seqan3::search_cfg::mode` has been replaced by individual configuration elements
    ([\#1639](https://github.com/seqan/seqan3/pull/1639)):
    * `seqan3::search_cfg::mode{seqan3::search_cfg::all}` to `seqan3::search_cfg::hit_all{}`
    * `seqan3::search_cfg::mode{seqan3::search_cfg::best}` to `seqan3::search_cfg::hit_single_best{}`
    * `seqan3::search_cfg::mode{seqan3::search_cfg::all_best}` to `seqan3::search_cfg::hit_all_best{}`
    * `seqan3::search_cfg::mode{seqan3::search_cfg::strata{5}}` to `seqan3::search_cfg::hit_strata{5}`
    * The `seqan3::search_cfg::hit_strata` member variable `value` has been replaced to `stratum`
  * The configuration`seqan3::search_cfg::output` has been replaced by individual configuration elements
    ([\#1862](https://github.com/seqan/seqan3/pull/1862)):
    * `seqan3::search_cfg::output{seqan3::search_cfg::text_position}` to `seqan3::search_cfg::output_reference_begin_position{}`
    * `seqan3::search_cfg::output{seqan3::search_cfg::text_position}` to `seqan3::search_cfg::output_index_cursor`
    * `seqan3::search_cfg::output_query_id{}` has been added
    * `seqan3::search_cfg::output_reference_id{}` has been added
* Removed `seqan3::bi_fm_index_cursor::to_rev_cursor()` and `seqan3::bi_fm_index::rev_cursor()`
  ([\#1892](https://github.com/seqan/seqan3/pull/1892)).

Header Changes:

```cpp
#include <seqan3/search/{algorithm => }/search.hpp>
#include <seqan3/{search/search_result_range.hpp => core/algorithm/algorithm_result_generator_range.hpp}>

#include <seqan3/search/configuration/{mode.hpp => hit.hpp}>
#include <seqan3/search/configuration/{max_error_rate.hpp => max_error.hpp AND max_error_common.hpp}>
```

## Notable Bug-fixes

#### Alignment

* When invoking the alignment algorithm with a user defined thread count using the `seqan3::align_cfg::parallel`
  configuration element, all available threads were used. This is now fixed and
  only the specified number of threads will be spawned ([\#1854](https://github.com/seqan/seqan3/pull/1854)).
* Using an unsigned score type via the `seqan3::align_cfg::score_type` configuration is prevented with a static assert,
  since gaps and mismatches have negative scores and thus need a signed score type
  ([\#1891](https://github.com/seqan/seqan3/pull/1891)).

#### Argument Parser

* Long option identifiers and their value must be separated by a space or equal sign `=`.
  Applying this restriction resolves an ambiguity that occurs if one long option identifier is the prefix of
  another ([\#1792](https://github.com/seqan/seqan3/pull/1792)).

  Valid short id value pairs: `-iValue`, `-i=Value`, `-i Value`
  Valid long id value pairs: `--id=Value`, `--id Value` (prohibited now: `--idValue`)

#### I/O

* The `seqan3::field::cigar` was added to the default fields for reading and writing alignment files
  ([\#1642](https://github.com/seqan/seqan3/pull/1642)).
  This has the following impact:
    1. Reading and writing in one line is now possible without additional reference information:
      `seqan3::alignment_file_output{"foo.sam"} = seqan3::alignment_file_input{"bar.sam"};`
    2. The `seqan3::alignment_file_output` now accepts `seqan3::field::cigar` and `seqan3::field::alignment`
       although they store redundant information. For the SAM/BAM format this ambiguity is handled by favouring the
       CIGAR information at all times if present.
  Note that this breaks your code if you have not selected custom fields and used structural bindings!

#### Search

* The `seqan3::fm_index_cursor::extend_right()`, `seqan3::bi_fm_index_cursor::extend_right()` and
  `seqan3::bi_fm_index_cursor::extend_left()` functions handle c-style strings without including the null character
  ([\#1588](https://github.com/seqan/seqan3/pull/1588)).
* `seqan3::fm_index` and `seqan3::bi_fm_index` construct the index correctly if a collection with a single text is
  passed as input ([\#1892](https://github.com/seqan/seqan3/pull/1892)).

#### Range

* Added `size()` function to `seqan3::views::kmer_hash`
  ([\#1722](https://github.com/seqan/seqan3/pull/1722)).
* `operator[](difference_type const n)` of the iterator of the `seqan3::views::kmer_hash` is declared `const`
  and returns value `n` steps after the current position without jumping to that position
  ([\#1756](https://github.com/seqan/seqan3/pull/1756)).

# 3.0.1

Note that 3.1.0 will be the first API stable release and interfaces in this release might still change.

## New features

#### Alphabet

* Added `seqan3::semialphabet_any`, a semi-alphabet that type erases all other semi-alphabets of the same size
  ([\#981](https://github.com/seqan/seqan3/pull/981)).
* Added `seqan3::dna3bs`, an alphabet that mimics a bisulfite-treated dna4 sequence
  ([\#1191](https://github.com/seqan/seqan3/pull/1191)).

#### Alignment

* The score type used in the alignment score matrix and the result type is configurable through a template
  argument of the `seqan3::align_cfg::result` configuration
  ([\#1340](https://github.com/seqan/seqan3/pull/1340)).
* The function`seqan3::align_pairwise` can be parallelised using the`seqan3::align_cfg::parallel` configuration
  ([\#1379](https://github.com/seqan/seqan3/pull/1379),
  [\#1444](https://github.com/seqan/seqan3/pull/1444)).

#### Argument parser

* Simplified reading file extensions from formatted files with the`seqan3::input_file_validator` and
  `seqan3::output_file_validator`
  ([\#863](https://github.com/seqan/seqan3/pull/863)).
* The `seqan3::value_list_validator` is now constructible from a range or a parameter pack
  ([\#1298](https://github.com/seqan/seqan3/pull/1298)).
* Enable subcommand argument parsing, see [How-to](https://docs.seqan.de/seqan/3.0.1/subcommand_arg_parse.html)
  for an example
  ([\#1185](https://github.com/seqan/seqan3/pull/1185)).
* The `seqan3::argument_parser::add_option` (and add_positional_option) calls allow enum types when using the
  `seqan3::enumeration_names` customisation point
  ([\#1196](https://github.com/seqan/seqan3/pull/1196)).

#### Build system

* `find_package(SeqAn3)` is now case-insensitive and always populates `SEQAN3_*` variables in all upper-case
  ([\#1427](https://github.com/seqan/seqan3/pull/1427)).

#### Core

* Added `seqan3::lzcnt`, `seqan3::tzcnt`, and `seqan3::popcount` for bit manipulation
  ([\#1141](https://github.com/seqan/seqan3/pull/1141)).
* Added traits for "metaprogramming" with `seqan3::type_list` and type packs
  ([\#1204](https://github.com/seqan/seqan3/pull/1204),
  [\#1214](https://github.com/seqan/seqan3/pull/1214),
  [\#1273](https://github.com/seqan/seqan3/pull/1273)).
* Added SIMD functions `seqan3::upcast` and `seqan3::upcast_signed`
  ([\#1190](https://github.com/seqan/seqan3/pull/1190)).

#### I/O

* We increased our input performance using a faster iterator on the stream buffer
  ([\#1030](https://github.com/seqan/seqan3/pull/1030)).
* Support of padded alignments in the SAM/BAM format was added
  ([\#1173](https://github.com/seqan/seqan3/pull/1173)).
* Reading `seqan3::field::cigar` into a vector over `seqan3::cigar` is supported via
  `seqan3::alignment_file_input`
  ([\#1192](https://github.com/seqan/seqan3/pull/1192)).
* Writing `seqan3::field::cigar` into a vector over `seqan3::cigar` is supported via
  `seqan3::alignment_file_output`
  ([\#1192](https://github.com/seqan/seqan3/pull/1192)).
* Asynchronous input (background file reading) supported via `seqan3::view::async_input_buffer`
  ([\#1205](https://github.com/seqan/seqan3/pull/1205)).

#### Range

* Added `seqan3::views::kmer_hash`, a view that computes hash values of an alphabet sequence given a
  `seqan3::shape`
  ([\#946](https://github.com/seqan/seqan3/pull/946)).
* Added `seqan3::views::to`, a view that returns a container created from a range by copying all elements
  ([\#1033](https://github.com/seqan/seqan3/pull/1033)).
* Added `seqan3::dynamic_bitset`, a container that stores single bits and has a dynamic size
  ([\#1153](https://github.com/seqan/seqan3/pull/1153)).
* Added `seqan3::views::translate_join`, analogue to `seqan3::views::translate` but returns a flattened range
  ([\#1171](https://github.com/seqan/seqan3/pull/1171)).
* Added `seqan3::views::to_simd`, a view that transforms a range of ranges into chunks of `seqan3::simd` vectors
  ([\#1190](https://github.com/seqan/seqan3/pull/1190)).
* Added `seqan3::views::as_const`, a view that provides only `const &` to elements of the underlying range
  ([\#1410](https://github.com/seqan/seqan3/pull/1410)).
* Added `seqan3::views::move`, a view that turns lvalue-references into rvalue-references
  ([\#1410](https://github.com/seqan/seqan3/pull/1410)).
* Renamed `seqan3::views::all` to `seqan3::views::type_reduce`
  ([\#1410](https://github.com/seqan/seqan3/pull/1410)).

#### Search

* The memory footprint of FM-indices over text collections was reduced
  ([\#1363](https://github.com/seqan/seqan3/pull/1363)).

#### Std

* We provide a `std::to_chars` overload for floating point data types in our `seqan3/std/from_chars` header
  ([\#1160](https://github.com/seqan/seqan3/pull/1160)).

## API changes

* The required version of the ranges-v3 library has increased:
  We now support the versions >= 0.10.0 and < 0.11.0, increasing the previous requirement of >= 0.5.0 and < 0.6.0
  ([\#1471](https://github.com/seqan/seqan3/pull/1471)).
* Customising for third party types has changes slightly:
  You are only affected if you added types to `seqan3::custom::`.
  Please see [About Customisation](https://docs.seqan.de/seqan/3.0.1/about_customisation.html)
  ([\#1225](https://github.com/seqan/seqan3/pull/1225)).
* All our concepts are named in the `snake_case` style
  (e.g. `seqan3::WritableAlphabet` -> `seqan3::writable_alphabet`)! This change was motivated by the decision of the
  ISO C++ committee to also use snake case everywhere
  ([\#1235](https://github.com/seqan/seqan3/pull/1235)).

#### Alphabet

* The `seqan3::cigar` alphabet is not an `seqan3::alphabet` anymore but only a `seqan3::semialphabet`
  ([\#1285](https://github.com/seqan/seqan3/pull/1285)).

#### Argument parser

* The `seqan3::value_list_validator` is not constructible from a std::initialiser_list anymore
  (e.g. `seqan3::value_list_validator{{1,2,3}}` does not work, use `seqan3::value_list_validator{1,2,3}` instead)
  ([\#1298](https://github.com/seqan/seqan3/pull/1298)).
* Changed class signature of input/output file validators:
  Most user code will be unaffected; to fix possible compiler errors you need to add an empty template parameter list to
  the respective instances (e.g. change `input_file_validator` to `input_file_validator<>`)
  ([\#863](https://github.com/seqan/seqan3/pull/863)).
* The member type that denotes which arguments a `validator` can validate has been renamed from `value_type` to
  `option_value_type`
  ([\#1394](https://github.com/seqan/seqan3/pull/1394)).
* Some exception names were altered and some removed ([\#1467](https://github.com/seqan/seqan3/pull/1467)):
  * The exception seqan3::parser_invalid_argument was renamed to seqan3::argument_parser_error.
  * The exception seqan3::validation_failed was renamed to seqan3::validation_error.
  * The exception seqan3::parser_design_error was renamed to seqan3::design_error and also inherits from
    seqan3::argument_parser_error.
  * The exception seqan3::type_conversion_failed was deprecated, you can catch seqan3::user_input_error instead.
  * The exception seqan3::overflow_error_on_conversion was deprecated, you can catch seqan3::user_input_error instead.

#### Build system

* [find_package](https://cmake.org/cmake/help/latest/command/find_package.html#version-selection) accepts minimum
  versions (e.g. `find_package(SEQAN3 3.0.1)` requires at least SeqAn3 with a version of `>= 3.0.1` and `< 4.0.0`)
  ([\#1425](https://github.com/seqan/seqan3/pull/1425)).
* The variable `SEQAN3_VERSION_STRING` defined by `find_package(SEQAN3)` was renamed to `SEQAN3_VERSION`
  ([\#1425](https://github.com/seqan/seqan3/pull/1425)).

#### Core

* The `type_list` header has moved:
  If you included `<seqan3/core/type_list.hpp>` you need to change the path to `<seqan3/core/type_list/type_list.hpp>`
  ([\#1204](https://github.com/seqan/seqan3/pull/1204)).

#### I/O

* Removed the field-based in- and output interface for sequence and structure files through std::get and std::tie:
  Output can instead be achieved with `seqan3::views:zip()`, for input we will implement `unzip()` in the future
  ([\#1398](https://github.com/seqan/seqan3/pull/1398)
  [\#1412](https://github.com/seqan/seqan3/pull/1412)).
* The `seqan3::field::flag` of SAM/BAM input and output is now an enum instead of an integer, see `seqan3::sam_flag`
  ([\#1390](https://github.com/seqan/seqan3/pull/1390)).
* Uppercase `seqan3::field` names are deprecated. Use the lower case field names instead. You can easily find
  and replace all occurrences by the following regex: find `field::([A-Z_]+)` replace `field::\L$1`
  ([\#1421](https://github.com/seqan/seqan3/pull/1421)).

* Removed the char type from the input and output files:
  Most user code will be unaffected; however, if you have fully specified all templates of any of the input or output
  files in your code, you need to remove the template parameter to select the char type of the stream,
  e.g. change `seqan3::sequence_file_input<traits_t, fields_t, formats_t, char>` to
  `seqan3::sequence_file_input<traits_t, fields_t, formats_t>`. Before this change, setting the char type gave the
  impression that also streams over wide characters are supported which is not the case yet
  ([\#1400](https://github.com/seqan/seqan3/pull/1400)).

#### Range

* The `seqan3::concatenated_sequences::data()` function has been deprecated:
  Use `seqan3::concatenated_sequences::raw_data()` instead
  ([\#1208](https://github.com/seqan/seqan3/pull/1208)).
* `seqan3::to_char` must always return a built-in character type
  ([\#1285](https://github.com/seqan/seqan3/pull/1285)).
* `seqan3/range/view` has be renamed to `seqan3/range/views`
  ([\#1251](https://github.com/seqan/seqan3/pull/1251)).
* namespace `seqan3::view` has been renamed to `seqan3::views`
  ([\#1251](https://github.com/seqan/seqan3/pull/1251)).

#### Search

* Changed class signature of (bi_)fm_index:
  All code that relies on automatic template deduction will be unaffected. In case you specified the template parameters
  of a `seqan3::fm_index` or `seqan3::bi_fm_index` you will need to add the alphabet type as first parameter and pass a
  `seqan3::text_layout` instead of a `bool` to indicate the text layout (single, collection).
  For example, `fm_index<false> index{text}` where `text` is of type `dna4_vector` needs to be changed to
  `fm_index<dna4, text_layout::single> index{text}`
  ([\#1222](https://github.com/seqan/seqan3/pull/1222)).

* The `construct()` method of the (bi_)fm_index is now private:
  Use the constructor `seqan3::fm_index::fm_index(text_t && text)` or `seqan3::bi_fm_index::bi_fm_index(text_t && text)`
  instead
  ([\#1222](https://github.com/seqan/seqan3/pull/1222)).

* The `seqan3::fm_index::char_type` member was renamed to `seqan3::fm_index::alphabet_type`
  The same applies for the `seqan3::bi_fm_index`
  ([\#1433](https://github.com/seqan/seqan3/pull/1433)).

* The `seqan3::fm_index_cursor::index_char_type` member was renamed to
  `seqan3::fm_index_cursor::index_alphabet_type`
  The same applies for the `seqan3::bi_fm_index_cursor`
  ([\#1433](https://github.com/seqan/seqan3/pull/1433)).

## Notable Bug-fixes

* All our headers are self contained
  ([\#1085](https://github.com/seqan/seqan3/pull/1085)).
* The alignment algorithm with edit distance returns the correct back coordinate
  ([\#1093](https://github.com/seqan/seqan3/pull/1093)).
* Inserting or deleting gaps into an empty `seqan3::gap_decorator` does not cause assert anymore
  ([\#1109](https://github.com/seqan/seqan3/pull/1109)).
* Some fixes to edge cases in BAM file writing
  ([\#1110](https://github.com/seqan/seqan3/pull/1110)).
* The application name of the `seqan3::argument_parser` is restricted to alpha-numeric characters and `_` and `-`
  ([\#1133](https://github.com/seqan/seqan3/pull/1133)).
* Copying and moving the `seqan3::fm_index` and `seqan3::bi_fm_index` now work properly
  ([\#1144](https://github.com/seqan/seqan3/pull/1144)).
* Searching in the `seqan3::fm_index` and `seqan3::bi_fm_index` constructed from a text collection containing a
  single text now returns the correct result
  ([\#1316](https://github.com/seqan/seqan3/pull/1316)).
* The view `seqan3::views::take` is sized if the underlying range is sized
  ([\#1146](https://github.com/seqan/seqan3/pull/1146)).
* The detection of the pthread library works correctly on linux based systems
  ([\#1200](https://github.com/seqan/seqan3/pull/1200)).
* The translation table for nucleotide to amino acid translation was corrected
  ([\#1485](https://github.com/seqan/seqan3/pull/1485)).
* The amino acid score matrices were corrected
  ([\#1455](https://github.com/seqan/seqan3/pull/1455)).

# 3.0.0 ("Escala")

This is the initial release of SeqAn3.
It is an entirely new library so there is no changelog that covers the differences to SeqAn2.

Note that 3.1.0 will be the first API stable release and interfaces in this release might still change.
