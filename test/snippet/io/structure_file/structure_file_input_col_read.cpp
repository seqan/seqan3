#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/alphabet/structure/wuss.hpp>
#include <seqan3/io/structure_file/input.hpp>
#include <seqan3/range/container/concatenated_sequences.hpp>

// Create a data struct used to hold the sequence information.
struct data_storage_t
{
    seqan3::concatenated_sequences<seqan3::rna5_vector>         sequences;
    seqan3::concatenated_sequences<std::string>                 ids;
    seqan3::concatenated_sequences<std::vector<seqan3::wuss51>> structures;
};

auto input = R"(> S.cerevisiae_tRNA-PHE M10740/1-73
GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA
(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)
> example
UUGGAGUACACAACCUGUACACUCUUUC
..(((((..(((...)))..)))))... (-3.71))";

int main()
{
    using seqan3::get;

    data_storage_t data_storage{};

    seqan3::structure_file_input fin{std::istringstream{input}, seqan3::format_vienna{}};

    // we move the buffer directly into our storage:
    data_storage.sequences  = std::move(get<seqan3::field::SEQ>(fin));
    data_storage.ids        = std::move(get<seqan3::field::ID>(fin));
    data_storage.structures = std::move(get<seqan3::field::STRUCTURE>(fin));
}
