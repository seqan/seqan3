# How to use cereal {#howto_use_cereal}

[TOC]

This HowTo documents how to use cereal to load from an archive or store into one. Every SeqAn data structure which is
marked as cerealisable can be used.

\tutorial_head{Easy, 15 min, , }

# Motivation

Storing and loading already existing indices is a common use case, thanks to the cereal library doing so is incredible
easy. This page will show you how to use cereal in SeqAn for one example. As example data structure std::vector was
picked, but as already mentioned any SeqAn data structure that is documented as cerealisable can be used.

# Storing

Storing a data structure is as easy as using the `cereal::BinaryOutputArchive`.

```cpp
#include <fstream>
#include <vector>

#if SEQAN3_WITH_CEREAL
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#endif

#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/std/filesystem>

// Written for std::vector, other types also work.
void store(std::vector<int16_t> const & data, std::filesystem::path tmp_dir)
{
    std::ofstream os(tmp_dir, std::ios::binary);      // Where output should be stored.
    cereal::BinaryOutputArchive archive(os);
    archive(data);                                    // Store data.
}

int main()
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path()/"data.out"; // Get the temp directory.
    store({1,2,3,4}, tmp_dir);                                                         // Calls store on a std::vector.

    std::filesystem::remove(tmp_file);                                                 // Remove the temporary file.
    return 0;
}
```

# Loading

Loading a data structure is as easy as using the `cereal::BinaryInputArchive`.

```cpp
#include <fstream>
#include <vector>

#if SEQAN3_WITH_CEREAL
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#endif

#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/std/filesystem>

// Written for std::vector, other types also work.
void load(std::vector<int16_t> & data, std::filesystem::path tmp_dir)
{
    std::ifstream is(tmp_dir, std::ios::binary);                      // Where input can be found.
    cereal::BinaryInputArchive archive(is);
    archive(data);                                                    // Load data.
}

// Written for std::vector, other types also work.
void store(std::vector<int16_t> const & data, std::filesystem::path tmp_dir)
{
    std::ofstream os(tmp_dir, std::ios::binary);      // Where output should be stored.
    cereal::BinaryOutputArchive archive(os);
    archive(data);                                    // Store data.
}

int main()
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path()/"data.out"; // Get the temp directory.

    store({1,2,3,4}, tmp_dir);                                                         // Calls store on a std::vector.
    std::vector<int16_t> vec;
    load(vec, tmp_dir);                                                                // Calls load on a std::vector.

    seqan3::debug_stream << vec << '\n';                                               // Prints [1,2,3,4].

    std::filesystem::remove(tmp_file);                                                 // Remove the temporary file.
    return 0;
}

```

# Storing & Loading in the same function

In the example above loading and storing was encapsulated in separated functions, it is possible to use
`cereal::BinaryInputArchive` and `cereal::BinaryOutputArchive` in one function, but it is necessary to encapsulate them then
with {}, otherwise one might encounter errors.
