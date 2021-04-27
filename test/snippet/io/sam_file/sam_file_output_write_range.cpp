#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sam_file/output.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::sam_file_output fout{std::ostringstream{}, seqan3::format_sam{}};

    std::vector<std::tuple<seqan3::dna5_vector, std::string>> range
    {
        { "ACGT"_dna5, "First" },
        { "NATA"_dna5, "2nd" },
        { "GATA"_dna5, "Third" }
    }; // a range of "records"

    fout = range; // will iterate over the records and write them
    // equivalent to:
    range | fout;
}
