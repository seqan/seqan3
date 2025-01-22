# First steps with SeqAn {#tutorial_first_example}

<!-- SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
     SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
     SPDX-License-Identifier: CC-BY-4.0
-->

***Learning Objective:***

This tutorial walks you through small SeqAn programs. It is intended to give you a short overview
of what to expect in the other tutorials and how to use this documentation.

\tutorial_head{Easy, 30 min, \ref setup, }

*Every page in the tutorials begins with the above table. It is recommended that you do the "prerequisite tutorials"
before the current one. You should also have a look at the links provided in "recommended reading" and maybe keep
them open in separate tabs/windows as reference.*

***These tutorials try to briefly introduce C++ features not well known. However, they do not teach programming in C++!
If you know how to program in another language, but are not familiar with C++ and/or the significant
changes in the language in recent years, we recommend the following resources:***

  * Bjarne Stroustrup: "A Tour of C++", Second Edition, 2018.

[TOC]

# Hello World!

Most good tutorials start with an easy *Hello World!* program. So have a look:

\snippet introduction_hello_world.cpp hello

\note
This is a code snippet. You will see many code snippets in our documentation.
Most of them are compilable as-is, but some are only valid in their context,
e.g. they depend on other code snippets given before/after the current one or
other statements implied by the text. You can **copy'n'paste** freely from these examples,
this implies no copyright-obligations (however distributing SeqAn or an application
using it does, see [Copyright](https://docs.seqan.de/seqan3/main_user/about_copyright.html) and
[Citing](https://docs.seqan.de/seqan3/main_user/about_citing.html)).

You may ask why we do not use std::cout or std::cerr for console output.
Actually, for the given text it does not make a difference since seqan3::debug_stream prints to std::cerr as well.
However, the debug stream provides convenient output for SeqAn's types as well as widely used data structures
(e.g. std::vector), which is especially helpful when you debug or develop your program
(that's where the name originates from).

\assignment{Assignment 1: Debug stream}
Write a program that creates a std::vector of type `int` and initialise the vector with a few values.
Then print the vector with seqan3::debug_stream. Does your program also work with std::cerr?
\endassignment
\solution
\snippet introduction_debug_stream.cpp debug
\endsolution

The above is an assignment with solution. You will find assignments in the tutorials to practise the discussed contents.
We believe that programming them will help you to memorise better and that it makes the tutorials more interesting and
interactive. The solutions provide the intended use; but often there are multiple ways to solve an assignment,
so don't worry too much if your solution is different from ours.

# API documentation

While the tutorials provide you with a walkthrough of some of our modules, the
[API documentation](https://docs.seqan.de/seqan3/main_user/modules.html) will be the go-to reference when you start
developing code with SeqAn.

Some helpful tips when browsing our documentation:

* You can search for seqan3 entities with the **search bar** in the top-right corner.
  E.g., start typing `debug_str` and the pop-up will suggest the `debug_stream` for you.
* If you don't have a specific entity you are searching for, the **landing pages** of each module are always a good
  start. E.g., the [Alphabet landing page](https://docs.seqan.de/seqan3/main_user/group__alphabet.html) first lists
  all submodules (Adaptation, Aminoacid, ...) and general alphabet-related seqan3 entities, followed by a detailed
  description of our alphabet module. Searching for keywords on this page might point you in the right direction.
* If you know you've seen some code snippet somewhere but don't remember where, have a look at our
  [cookbook](https://docs.seqan.de/seqan3/main_user/cookbook.html). It is not structured and huge, but works
  well if you do a key word search with `Ctrl+F`.

We recommend you to open the API documentation in separate browser tab s.t. you can easily switch back to the tutorial.

If you have troubles or the documentation is missing some information, feel free to write to the developers
of SeqAn on [Github](https://github.com/seqan/seqan3/issues/new/choose) and ask your questions directly.

# Modules in SeqAn

Modules structure the SeqAn library into logical units. There are, for instance,

* [`alphabet`](https://docs.seqan.de/seqan3/main_user/group__alphabet.html): `seqan3::dna4` etc.
* [`io`](http://docs.seqan.de/seqan3/main_user/group__io.html): read/write FASTA, SAM, ...
* [`alignment`](http://docs.seqan.de/seqan3/main_user/group__alignment.html): compute pairwise alignments etc.
* [`search`](http://docs.seqan.de/seqan3/main_user/group__search.html): search via an FM-Index etc.

and some more.

Some modules consist of submodules and the module structure is represented by the file hierarchy in the `include/`
directory. Whenever you use functions of a module, make sure to `#include` the correct header file.

Each directory in the SeqAn sources contains an `all.hpp` file which includes all the functionality
of the respective module.
For small examples and quick prototyping, you can just include these `all.hpp`-headers.
However, for larger projects, we recommend you include only the necessary headers, because this will reduce the
compile time measurably.

\assignment{Assignment 2: Modules and API documentation}
In your program of assignment 1, initialise a vector of `seqan3::dna4` instead of `int`.
The vector shall store the DNA string `ACTG`.
Check the [API documentation](http://docs.seqan.de/seqan3/main_user/modules.html) for which header you need to include.
Additionally, browse the documentation for `seqan3::dna4` on how to initialise a `seqan3::dna4` letter.
\endassignment
\solution
\snippet introduction_dna4.cpp debug
\endsolution

# Some general notes that might help to dive into SeqAn

## SeqAn and the STL

In contrast to the former version of SeqAn (2.x.x releases), we try to be very close to the standard and all other
data structures and algorithms should work on STL data structures as well.

Analogous to the STL, SeqAn3 uses `snake_case` everywhere.

## Modern C++

We use a lot of Modern C++ in SeqAn, so some things might look alien at first,
e.g. type templates are used like ordinary types in many situations (no `<>`).
We also always use `{}` to initialise objects and not `()` which is only used for function calls.
In general, the style should be much easier for newcomers.

## Avoid using namespace seqan3

In concordance with the [C++ Core guidelines](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#Rs-using),
we encourage you to avoid declaring `using namespace seqan3;`. This has the benefit of easily distinguishing
between `seqan3` features and standard C++ (`std`). The only exception are string literals, where we often use
`using namespace seqan3::literals;` for convenience.

# The next tutorials

Now that you reached the end of this first tutorial, you know how SeqAn code looks like and you are able
to write some first code fragments. Let's go more into detail with the module-based tutorials!
