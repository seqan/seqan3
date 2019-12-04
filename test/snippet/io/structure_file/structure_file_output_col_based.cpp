#include <sstream>
#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/alphabet/structure/wuss.hpp>
#include <seqan3/io/structure_file/output.hpp>
#include <seqan3/range/container/concatenated_sequences.hpp>
#include <seqan3/range/views/zip.hpp>

using seqan3::operator""_rna5;
using seqan3::operator""_wuss51;

struct data_storage_t
{
    seqan3::concatenated_sequences<seqan3::rna5_vector>         sequences{"AACGUU"_rna5};
    seqan3::concatenated_sequences<std::string>                 ids{std::string{"seq1"}};
    seqan3::concatenated_sequences<std::vector<seqan3::wuss51>> structures{".(())."_wuss51};
};

int main()
{
    data_storage_t data_storage{}; // a global or globally used variable in your program

    // ... in your file writing function:

    seqan3::structure_file_output fout{std::ostringstream{}, seqan3::format_vienna{}};

    fout = seqan3::views::zip(data_storage.sequences, data_storage.ids, data_storage.structures);
}
