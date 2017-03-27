# SeqAn3 -- the modern C++ library for sequence analysis

SeqAn3 is the next version of the popular SeqAn template library for the analysis of biological sequences. If you develop applications, we recommend you currently stick to [SeqAn2](https://github.com/seqan/seqan).

SeqAn3 is major redesign of SeqAn and has fundamental changes throughout the codebase. We expect that porting applications from SeqAn2 to SeqAn3 requires substantial work, however by embracing new technologies from C++17 and modern third party libraries, SeqAn3 will be a much improved experience over SeqAn2.


## Quick facts

  * same design goals as SeqAn2: fast, efficient, extensible C++ header library for sequence analysis
  * different design patterns: generic programming via C++ Concepts, encapsulation and members
  * modern C++ that relies heavily on C++17, the [Concepts TS](http://www.stroustrup.com/good_concepts.pdf) and [Ranges TS](https://github.com/ericniebler/range-v3)


## Requirements

### Users of the library

To include SeqAn3 in your app, you need the following:

|                   | requirement                                          | version  | comment                                   |
|-------------------|------------------------------------------------------|----------|-------------------------------------------|
|**compiler**       | [GCC](http://gcc.gnu.org)                            | ≥ 7      | no other compiler is currently supported! |
|**build system**   | [cmake](https://cmake.org)                           | ≥ 3.4    | optional, but recommended                  |
|**required libs**  | [SDSL](https://github.com/xxsds/sdsl-lite)           | ≥ 3      | succint datastructures                    |
|                   | [Ranges-V3](https://github.com/ericniebler/range-v3) | == 0.2.1 | ranges and views                          |
|**optional libs**  | [Cereal](https://github.com/USCiLab/cereal)          | ≥ 1.2    | serialization                              |
   
### Developers of the library

To build the tests and API documentation, you also need:

|                   | requirement                                          | version  |
|-------------------|------------------------------------------------------|----------|
|**build system**   | [cmake](https://cmake.org)                           | ≥ 3.4    |
|**test system**    | [GoogleTest](https://github.com/google/googletest)   | ≥ 1.8    | 
|**doc system**     | [Doxygen](https://github.com/doxygen/doxygen)        | ≥ 1.8    |

## Using

To use SeqAn3, clone this repository including it's submodules:

```sh
git clone --recursive https://github.com/seqan/seqan3.git
```

SeqAn3 is still a header-only library so it is sufficient to add the `include` folder to your include path. When building make sure that you add the required parameters:

```sh
g++-7 -std=c++17 -fconcepts -I /path/to/seqan3/include -I /path/to/seqan3/range-v3/include -I /path/to/seqan3/sdsl-lite/include myfile.cpp
```

## Documentation

### Recommended reading for (modern) C++

  * [C++ Core Guidelines](https://github.com/isocpp/CppCoreGuidelines/blob/master/CppCoreGuidelines.md)
  
### Users of the library

  * [API documentation](https://seqan3-api.readthedocs.org)
  * [Manual and tutorials](https://seqan3-manual.readthedocs.org)

  
### Developers of the library

  * [Github Wiki](https://github.com/seqan/seqan3/wiki)
  
