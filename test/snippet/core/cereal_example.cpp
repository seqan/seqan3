#include <vector>
#include <fstream>

#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/std/filesystem>

void load(std::vector<int> & data, std::filesystem::path tmp_dir) // Written for std::vector, other types also work
{
    std::ifstream is(tmp_dir/"data.out", std::ios::binary);       // Where input can be found
    cereal::BinaryInputArchive archive(is);
    archive(data);                                                // load data
}

void store(std::vector<int> data, std::filesystem::path tmp_dir) // Written for std::vector, other types also work
{
    std::ofstream os(tmp_dir/"data.out", std::ios::binary);      // Where output should be stored
    cereal::BinaryOutputArchive archive(os);
    archive(data);                                               // Store data
}

int main()
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory

    // The following example is for an std::vector but any seqan3 data structure that is documented as serialisable
    // could be used, e.g. fm_index.
    store({1,2,3,4}, tmp_dir);                                              // Calls store on a std::vector,
    std::vector<int> vec;
    load(vec, tmp_dir);                                                     // Calls load on a std::vector

    seqan3::debug_stream << vec << "\n";                                    // Prints [1,2,3,4]
    return 0;
}
