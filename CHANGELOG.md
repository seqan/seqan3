# Changelog {#about_changelog}

[TOC]

This changelog contains a top-level entry for each release with sections on new features, API changes and notable
bug-fixes (not all bug-fixes will be listed).

Get to know SeqAn3 with our [tutorials](https://docs.seqan.de/seqan/3-master-user/usergroup1.html).

Please see the release announcement: https://www.seqan.de/announcing-seqan3/

See the porting guide for some help on porting: https://docs.seqan.de/seqan/3-master-user/howto_porting.html

See the documentation on [API stability](https://docs.seqan.de/seqan/3-master-user/about_api.html) to learn about
when API changes are allowed.

<!--
The following API changes should be documented as such:
  * a previously experimental interface now being marked as stable
  * an interface being removed
  * syntactical changes to an interface (e.g. renaming or reordering of files, functions, parameters)
  * semantic changes to an interface (e.g. a function's result is now always one larger) [DANGEROUS!]

If possible, provide tooling that performs the changes, e.g. a shell-script.
-->

# 3.0.3

## New features

#### Alphabet

* Added `seqan3::phred94`, a quality type that represents the full Phred Score range (Sanger format) and is used for
  PacBio Phred scores of HiFi reads ([\#2290](https://github.com/seqan/seqan3/pull/2290)).

#### Argument Parser

* We expanded the `seqan3::output_file_validator`, with a parameter `seqan3::output_file_open_options` to allow overwriting
  output files ([\#2009](https://github.com/seqan/seqan3/pull/2009)).
* The `seqan3::argument_parser` has a new member function `seqan3::argument_parser::is_option_set` that
  checks whether an option, identified by its long or short name, was set on the command line by the user
  ([\#1859](https://github.com/seqan/seqan3/pull/1859)).
* The subcommand of the `seqan3::argument_parser` may contain a dash
  ([\#2306](https://github.com/seqan/seqan3/pull/2306)).

#### Build system

* We now use Doxygen version 1.9.1 to build our documentation ([\#2327](https://github.com/seqan/seqan3/pull/2327)).

#### Search

* The `seqan3::fm_index_cursor` exposes its suffix array interval ([\#2076](https://github.com/seqan/seqan3/pull/2076)).

## Notable Bug-fixes

### Argument Parser

* Argument parsing with subcommands: If the user forgets or misspells the subcommand, the error is thrown when calling
  `seqan3::argument_parser::parse()` and not on construction of the `seqan3::argument_parser`
  ([\#2179](https://github.com/seqan/seqan3/pull/2179)).
* The `seqan3::regex_validator` parses `std::filesystem::path`'s correctly now
  ([\#2216](https://github.com/seqan/seqan3/pull/2216)).

## API changes

#### Alphabet

* Removed seqan3::char_is_valid_for requirement from seqan3::writable_alphabet and seqan3::writable_constexpr_alphabet
  ([\#2290](https://github.com/seqan/seqan3/pull/2337)).

#### Argument Parser

* The enum names of `seqan3::option_spec` were changed to lower case
  ([\#2285](https://github.com/seqan/seqan3/pull/2285)):
  * `seqan3::option_spec::DEFAULT` is replaced by `seqan3::option_spec::standard`.
  * `seqan3::option_spec::REQUIRED` is replaced by `seqan3::option_spec::required`.
  * `seqan3::option_spec::ADVANCED` is replaced by `seqan3::option_spec::advanced`.
  * `seqan3::option_spec::HIDDEN` is replaced by `seqan3::option_spec::hidden`.

# 3.0.2

Note that 3.1.0 will be the first API stable release and interfaces in this release might still change.

* Check out our new [SeqAn3 Cookbook](https://docs.seqan.de/seqan/3.0.2/cookbook.html). It contains a listing of code
  examples on how to perform particular tasks using the library.

* SeqAn 3.0.2 is known to compile with GCC 7.5, 8.4, 9.3 and 10.2. Future versions (e.g. GCC 10.3 and 11) might work,
  but weren’t yet available at the time of this release.

## Notable Bug-fixes

#### I/O

* The SAM/BAM format reads the quality field (or any other text field) that starts with an asterisk (`*`) but is not
  empty correctly now ([\#2184](https://github.com/seqan/seqan3/pull/2184)).

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
  `seqan3::semialphabet` ([\#1721](https://https://github.com/seqan/seqan3/pull/1721)).

#### Search

* Added `seqan3::interleaved_bloom_filter`, a data structure that efficiently answers set-membership queries for
  multiple bins ([\#920](https://github.com/seqan/seqan3/pull/920)).
* Added `seqan3::search_cfg::hit`, which allows dynamic configuration of the hit strategy.
  ([\#1853](https://github.com/seqan/seqan3/pull/1853)).
* Added `seqan3::search_cfg::on_result`, which allows providing a custom callback for the search algorithm
  ([\#2019](https://github.com/seqan/seqan3/pull/2019)).
* The `seqan3::fm_index_cursor` and `seqan3::bi_fm_index_cursor` can be serialised
  ([\#2048](https://github.com/seqan/seqan3/pull/2019)).

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
  ([\#2025](https://https://github.com/seqan/seqan3/pull/2025)).

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

* The `seqan3::value_list_validator` is not constructible from a std::initialiser_list any more
  (e.g. `seqan3::value_list_validator{{1,2,3}}` does not work, use `seqan3::value_list_validator{1,2,3}` instead)
  ([\#1298](https://github.com/seqan/seqan3/pull/1298)).
* Changed class signature of input/output file validators:
  Most user code will be unaffected; to fix possible compiler errors you need to add an empty template parameter list to
  the respective instances (e.g. change `input_file_validator` to `input_file_validator<>`)
  ([\#863](https://github.com/seqan/seqan3/pull/863)).
* The member type that denotes which arguments a `validator` can validate has been renamed from `value_type` to
  `option_value_type`
  ([\#1394](https://github.com/seqan/seqan3/pull/1394)).
* Some exception names were altered and some removed ([\#1394](https://github.com/seqan/seqan3/pull/1467)):
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
