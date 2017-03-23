# SeqAn3 -- the modern C++ library for sequence analysis

SeqAn3 is the next version of the popular SeqAn template library for the analysis of biological sequences. If you develop applications, we recommend you currently stick to [SeqAn2](https://github.com/seqan/seqan).

SeqAn3 is major redesign of SeqAn and has fundamental changes throughout the codebase. We expect that porting applications from SeqAn2 to SeqAn3 requires substantial work, however by embracing new technologies from C++17 and modern third party libraries, SeqAn3 will be a much improved experience over SeqAn2.


## Quick facts

* modern C++17 library that uses C++ Concepts and C++ Ranges
* **requires GCC >= 7**, other compilers are currently not supported!
* API that is closer to traditional object-oriented C++ with members and less (visible) templates
* makes use of the following modern C++ library:
   * [SDSL](https://github.com/xxsds/sdsl-lite)
   * [Ranges-V3](https://github.com/ericniebler/range-v3)
   * [Cereal](https://github.com/USCiLab/cereal)
   

## Using

To use SeqAn3, clone this repository including it's submodules:

```sh
git clone --recursive https://github.com/seqan/seqan3.git
```

SeqAn3 is still a header-only library so it is sufficient to add the `include` folder to your include path. When building make sure that you add the required parameters:

```sh
g++-7 -std=c++17 -fconcepts -I /path/to/seqan3/include -I /path/to/seqan3/range-v3/include -I /path/to/seqan3/sdsl-lite/include myfile.cpp
```

