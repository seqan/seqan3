// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <fstream>
#include <vector>

//! [binary_include]
#include <cereal/archives/binary.hpp> // includes the cereal::BinaryOutputArchive
//! [binary_include]
//! [vector_include]
#include <cereal/types/vector.hpp> // includes cerealisation support for std::vector
//! [vector_include]

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/test/tmp_directory.hpp>

// Written for std::vector, other types also work.
void store(std::vector<int16_t> const & data, std::filesystem::path const & tmp_file)
{
    std::ofstream os(tmp_file, std::ios::binary); // Where output should be stored.
    cereal::BinaryOutputArchive archive(os);      // Create an output archive from the output stream.
    archive(data);                                // Store data.
}

int main()
{
    // The following example is for a std::vector but any seqan3 data structure that is documented as serialisable
    // could be used, e.g. seqan3::fm_index.
    seqan3::test::tmp_directory tmp{};
    auto tmp_file = tmp.path() / "data.out"; // This is a temporary file name, use any other filename.

    std::vector<int16_t> vec{1, 2, 3, 4};
    store(vec, tmp_file); // Calls store on a std::vector.

    return 0;
}
