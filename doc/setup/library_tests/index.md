# Setting up library tests {#setup_tests}

<!-- SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
     SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
     SPDX-License-Identifier: CC-BY-4.0
-->

<b>Learning Objective:</b><br>
In this guide you will learn to set up SeqAn3's library tests to make sure that your contributions don't break anything.

\tutorial_head{Easy, 20 Minutes, \ref setup, }

[TOC]

# Unit testing

Unit tests are the most important tests, they should cover all functionality and all code-paths of the library.
Before submitting a pull request to our repository, make sure that the unit tests pass on your system.

## Setting up unit tests

Assume that you have cloned SeqAn into `/home/me/devel/seqan3` and performed some local changes.

Create an out-of-source build directory and change to it:

```bash
mkdir -p /home/me/devel/seqan3-build/debug
cd /home/me/devel/seqan3-build/debug
```

Invoke CMake (this is often referred to as the "Configure" step):

```bash
cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_COMPILER=g++-11 ../../seqan3/test/unit
```

The build type could be "Release", but stick to "Debug" if you aim to make contributions to the SeqAn3 codebase.
Specifying the compiler is optional; depending on your setup, you may need to give the full path.

### Building all unit tests

Invoke Make (this is the actual build step):

```bash
make
```

You may want to pass `-j x` where is `x` is the number of CPUs on your system to speed up the build.
If you have build failures, it is recommend to run make again without `-j x` to receive better readable output.

After all test have been successfully built, run ctest to ensure correct results:

```bash
ctest .
```

Fix any issues you encounter and re-run make and ctest.
You only need to re-run cmake if you add or remove unit tests in between your changes.

### Building a specific unit test

If you are working on a very specific piece of code, it might be beneficial to ensure that specific test builds and
passes before re-building the entire set.

This builds the test for seqan3::dna4:

```bash
make dna4_test
```

And this runs the test:

```bash
alphabet/nucleotide/dna4_test
```

Note that you need to give the (relative) path when running the executable, but not when building the test.
Running the test executable individually will also tell you which parts of the test fail.

\attention Before you (re-)submit changes in a pull-request, please build and run **all unit tests** as small changes
can have unexpected side-effects.

# Other test suites

SeqAn has the following test suites:

  * unit: tests the API of the library
  * documentation: tests that everything is properly documented
  * snippet: tests the buildability of code snippets inside the documentation
  * header: tests that every header includes all required headers and detects linkage issues
  * performance: contains microbenchmarks
  * macro_benchmark: contains macrobenchmarks

If your unit tests pass and you submit a pull request, our continuous integration builds also run the other test suites.
You might encounter failures in one of them in which case you need to also setup that test suite to fix your problems
before you re-submit/update your pull request.

The only difference in setting up the other test suites is a different path when invoking CMake.
Since you used out-of-source builds, you can simply create another directory for the other test suite(s).

This will setup snippet tests:

```bash
mkdir -p /home/me/devel/seqan3-build/snippet
cd /home/me/devel/seqan3-build/snippet

cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_COMPILER=g++-11 ../../seqan3/test/snippet

make -j 4

ctest .
```

Documentation tests do not require setting a build type or compiler, but they require that doxygen be installed on
the system. Locally built documentation will be placed into
`/home/me/devel/seqan3-build/documentation/doc_usr/html/`.

# Platform specific notes

  * On *BSD operating systems (but not macOS), call `gmake` whenever you would call `make`.
