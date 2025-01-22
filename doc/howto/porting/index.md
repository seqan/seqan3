# Porting from SeqAn2 {#howto_porting}

<!-- SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
     SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
     SPDX-License-Identifier: CC-BY-4.0
-->

SeqAn3 is a completely new library so there is no easy or automated way of porting applications from SeqAn2.
We recommend that everyone go through the tutorial, learn the new idioms and then re-write their application.
We know this is a big burden and we promise that next major versions will be incremental updates again.
The notes on this page may help transitioning and/or maintaining dependencies to SeqAn2 and SeqAn3.

[TOC]

# Important differences to prior versions

Please see \ref tutorial_concepts and \ref tutorial_ranges for some of the fundamental, structural differences
and design choices of SeqAn3.

## Naming

SeqAn3 uses `snake_case` instead of `CamelCase` for naming all entities other than concepts.
Concepts are named using `CamelCase`. This should be very close to the standard library.

## Initialisation

SeqAn3 uses brace initialisation, also known as "uniform initialisation" for all class types.
Built-in arithmetic types like `uint64_t` may be initialised via `=`.

There are many online resources on "uniform initialisation", e.g.:

  * https://arne-mertz.de/2015/07/new-c-features-uniform-initialization-and-initializer_list/

## Return values and parameters

SeqAn1 and SeqAn2 used out-parameters as the primary way of returning values, SeqAn3 uses plain return values.
This is no less efficient nowadays, because there is guaranteed copy elision since C++17.

The details of copy elision are described here:

  * https://en.cppreference.com/w/cpp/language/copy_elision

# Using SeqAn3 and SeqAn2 in the same application

SeqAn2 and SeqAn3 can coexist in the same application, however certain care should be taken.

## Versions

There are some incompatibilities between SeqAn-2.4 and SeqAn3, please update to the main branch of SeqAn2 or
use the 2.5 release (when available).

## Paths and Namespaces

For SeqAn1 and SeqAn2 the following entities are called "seqan":

  * the primary namespace
  * the include folder
  * the CMake library name

All of the above are named "seqan3" for the SeqAn3 library so coexistence should be possible on all levels.

We discourage having any `using namespace seqan;` or `using namespace seqan3;` statements in your application and
and suggest replacing any occurrences of the former with fully qualified lookups, e.g. writing `seqan::readRecord`.
Clearly marking all SeqAn2-code as such will help identify (and replace) it later on and is a good step before
introducing SeqAn3 into your code.
We are working on a script to help with this.

## Mixing types

Many SeqAn2 interfaces handle standard library containers and strings and we are working on adding a support layer
to SeqAn2 that handles SeqAn3-alphabets.
For now you may have to convert data structures if using old and new interfaces within the same application.

Some components that are rather independent and good first targets to start porting are:

  * The argument parser.
  * Input and Output.

## The argument parser version check feature

The major change is that we do not provide a cmake-level directive to turn off version checking.
Instead, the current implementation lets you disable the version check in the following ways:

  1. The developer can disable the feature permanently on construction of the seqan3::argument_parser.
  2. The administrator or user can deactivate the feature globally on their system via the environment variable
     `SEQAN3_NO_VERSION_CHECK`.
  3. The user can disable the feature temporarily for one day via a command-line option when calling the application.

Otherwise, if non of the above happened, the user will be prompted at the first call of the application to select
a configuration (Never/Always/Yes/No). However, if the application is executed in a non-terminal process,
e.g within a cron job or a workflow, no prompt will be issued but instead the version check will be triggered in
the background. This will happen no more than once per day.

If you are interested in the new argument parser design, take a look at our tutorial \ref tutorial_argument_parser.

## SeqAn2 features missing from SeqAn3

This is a non-exhaustive list of features that are currently missing, but **will be ported**:

  * Annotation file formats (GFF, GTF, VCF...)
  * Blast tabular file formats
  * Memory Mapping support
  * Journaled String Tree
  * k-mer/q-gram based indexing/searching

This is a non-exhaustive list of features that we **do not plan to port:**

  * The Store module
  * Online-Search algorithms (not used in practice)
  * The Graph module (use the lemon library instead)
  * Uncompressed suffix arrays (no benefit over FM indexes)
