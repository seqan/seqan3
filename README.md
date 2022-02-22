# SeqAn3 -- the modern C++ library for sequence analysis

[![build status][1]][2]
[![codecov][3]][4]
[![license][5]][6]
[![latest release][7]][8]
[![platforms][9]][10]
[![start][11]][12]
[![twitter][13]][14]

<!--
    Above uses reference-style links with numbers.
    See also https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet#links.

    For example, `[![build status][1]][2]` evaluates to the following:
        `[link_text][2]`
        `[2]` is a reference to a link, i.e. `[link_text](https://...)`

        `[link_text]` = `[![build status][1]]`
        `[1]` is once again a reference to a link - this time an image, i.e. `[![build status](https://...)]
        `![build status]` is the text that should be displayed if the linked resource (`[1]`) is not available

    `[![build status][1]][2]` hence means:
    Show the picture linked under `[1]`. In case it cannot be displayed, show the text "build status" instead.
    The picture, or alternative text, should link to `[2]`.
-->

[1]: https://img.shields.io/github/workflow/status/seqan/seqan3/SeqAn3%20CI%20on%20Linux/master?style=flat&logo=github&label=SeqAn3%20CI "Open GitHub actions page"
[2]: https://github.com/seqan/seqan3/actions?query=branch%3Amaster
[3]: https://codecov.io/gh/seqan/seqan3/branch/master/graph/badge.svg?token=BH1FQiBBle "Open Codecov page"
[4]: https://codecov.io/gh/seqan/seqan3
[5]: https://img.shields.io/badge/license-BSD-green.svg "Open Copyright page"
[6]: https://docs.seqan.de/seqan/3-master-user/about_copyright.html
[7]: https://img.shields.io/github/release/seqan/seqan3.svg "Get the latest release"
[8]: https://github.com/seqan/seqan3/releases/latest
[9]: https://img.shields.io/badge/platform-linux%20%7C%20bsd%20%7C%20osx-informational.svg "Read more about our API"
[10]: https://docs.seqan.de/seqan/3-master-user/about_api.html
[11]: https://img.shields.io/github/stars/seqan/seqan3.svg?style=social "See who starred us"
[12]: https://github.com/seqan/seqan3/stargazers
[13]: https://img.shields.io/twitter/follow/SeqAnLib.svg?label=follow&style=social "Follow us on Twitter"
[14]: https://twitter.com/seqanlib

SeqAn3 is the new version of the popular SeqAn template library for the analysis of biological sequences.
It enables the rapid development of high-performance solutions by providing generic algorithms and data structures
for:

  * sequence representation and transformation
  * full-text indexing and efficient search
  * sequence alignment
  * input/output of common file formats

By leveraging *Modern C++* it provides unprecedented ease-of-use without sacrificing performance.

Please see the [online documentation](https://docs.seqan.de/seqan/3-master-user/) for more details.

## Quick facts

  * C++ header-only library: easy to integrate with your app & easy to distribute
  * liberal open source license: allows integration with any app or library, requires only attribution
  * very high code quality standards: >97% unit test coverage, performance regression tests, ...
  * extensive API documentation & tutorials: more lines of documentation than lines of code
  * aims to support any 64-bit architecture running Linux/POSIX; currently big-endian CPU architectures
    like s390x are less supported

## Dependencies

|                   | requirement                                          | version  | comment                                     |
|-------------------|------------------------------------------------------|----------|---------------------------------------------|
|**compiler**       | [GCC](https://gcc.gnu.org)                           | ≥ 7      | no other compiler is currently supported!   |
|**build system**   | [CMake](https://cmake.org)                           | ≥ 3.4    | optional, but recommended                   |
|**required libs**  | [SDSL](https://github.com/xxsds/sdsl-lite)           | ≥ 3      |                                             |
|                   | [Range-V3](https://github.com/ericniebler/range-v3)  | ≥ 0.11.0 |                                             |
|**optional libs**  | [cereal](https://github.com/USCiLab/cereal)          | ≥ 1.2.3  | required for serialisation and CTD support  |
|                   | [zlib](https://github.com/madler/zlib)               | ≥ 1.2    | required for `*.gz` and `.bam` file support |
|                   | [bzip2](https://www.sourceware.org/bzip2)            | ≥ 1.0    | required for `*.bz2` file support           |

## Usage

We recommend that you use CMake to build your project:

  * [Setup-Tutorial](https://docs.seqan.de/seqan/3-master-user/setup.html)
  * Using CMake guarantees that all optional dependencies are automatically detected and activated.

Quick-Setup without CMake:

  * Clone the repository with submodules: `git clone --recurse-submodules https://github.com/seqan/seqan3.git`
  * Add the following to your compiler invocation:
    * the include directories of SeqAn and its dependencies
    * C++17 mode with concepts support
    * Macros indicating the presence of zlib and bzip2 (set only if actually available in your paths!)
  * The command could look like this:
```sh
g++-7 -O3 -DNDEBUG -Wall -Wextra                                \
    -std=c++17 -fconcepts                                       \
    -I       /path/to/seqan3/include                            \
    -isystem /path/to/seqan3/submodules/range-v3/include        \
    -isystem /path/to/seqan3/submodules/sdsl-lite/include       \
    -isystem /path/to/seqan3/submodules/cereal/include          \
    -DSEQAN3_HAS_ZLIB=1 -DSEQAN3_HAS_BZIP2=1                    \
    -lz -lbz2 -lstdc++fs -pthread                               \
  your_file.cpp
```

## Sponsorships

[![Vercel](https://raw.githubusercontent.com/seqan/seqan3/master/test/documentation/.vercel/powered-by-vercel.svg)](https://vercel.com/?utm_source=seqan&utm_campaign=oss)

Vercel is kind enough to sponsor our documentation preview-builds within our pull requests. Check them out!
