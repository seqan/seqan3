#if SEQAN3_WITH_CEREAL
#    include <fstream>
#    include <vector>

#    include <seqan3/core/debug_stream.hpp>
#    include <seqan3/test/tmp_directory.hpp>

#    include <cereal/archives/binary.hpp> // includes the cereal::BinaryInputArchive and cereal::BinaryOutputArchive
#    include <cereal/types/vector.hpp>    // includes cerealisation support for std::vector

// Written for std::vector, other types also work.
void load(std::vector<int16_t> & data, std::filesystem::path const & tmp_file)
{
    std::ifstream is(tmp_file, std::ios::binary); // Where input can be found.
    cereal::BinaryInputArchive archive(is);       // Create an input archive from the input stream.
    archive(data);                                // Load data.
}

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
    // could be used, e.g. fm_index.
    seqan3::test::tmp_directory tmp{};
    auto tmp_file = tmp.path() / "data.out"; // this is a temporary file path, use any other filename.

    std::vector<int16_t> vec{1, 2, 3, 4};
    store(vec, tmp_file); // Calls store on a std::vector.
    // This vector is needed to load the information into it.
    std::vector<int16_t> vec2;
    load(vec2, tmp_file); // Calls load on a std::vector.

    seqan3::debug_stream << vec << '\n'; // Prints [1,2,3,4].

    return 0;
}
#endif
