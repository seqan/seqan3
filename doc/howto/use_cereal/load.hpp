#include <fstream>
#include <vector>

//! [binary_include]
#include <cereal/archives/binary.hpp> // includes the cereal::BinaryInputArchive and cereal::BinaryOutputArchive
//! [binary_include]
#include <cereal/types/vector.hpp>    // includes cerealisation support for std::vector

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/test/tmp_filename.hpp>

// Written for std::vector, other types also work.
void load(std::vector<int16_t> & data, seqan3::test::tmp_filename & tmp_file)
{
    std::ifstream is(tmp_file.get_path(), std::ios::binary); // Where input can be found.
    cereal::BinaryInputArchive archive(is);                  // Create an input archive from the input stream.
    archive(data);                                           // Load data.
}

// Written for std::vector, other types also work.
void store(std::vector<int16_t> const & data, seqan3::test::tmp_filename & tmp_file)
{
    std::ofstream os(tmp_file.get_path(), std::ios::binary); // Where output should be stored.
    cereal::BinaryOutputArchive archive(os);                 // Create an ouput archive from the output stream.
    archive(data);                                           // Store data.
}

int main()
{
    // The following example is for an std::vector but any seqan3 data structure that is documented as serialisable
    // could be used, e.g. seqan3::fm_index.
    seqan3::test::tmp_filename tmp_file{"data.out"}; // This is a temporary file, use any other filename.

    std::vector<int16_t> vec{1,2,3,4};
    store(vec, tmp_file);                            // Calls store on a std::vector.
    // This vector is needed to load the information into it.
    std::vector<int16_t> vec2;
    load(vec2, tmp_file);                            // Calls load on a std::vector.

    seqan3::debug_stream << vec2 << '\n';             // Prints [1,2,3,4].

    return 0;
}
