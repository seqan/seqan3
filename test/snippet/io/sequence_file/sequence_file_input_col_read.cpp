#include <sstream>
#include <string>
#include <utility>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/range/container/concatenated_sequences.hpp>
#include <seqan3/range/view/get.hpp>

auto input = R"(> TEST1
ACGT
> Test2
AGGCTGA
> Test3
GGAGTATAATATATATATATATAT)";

struct data_storage_t
{
    seqan3::concatenated_sequences<seqan3::dna5_vector>  sequences;
    seqan3::concatenated_sequences<std::string>          ids;
};

int main()
{
    using seqan3::get;

    data_storage_t data_storage{};

    seqan3::sequence_file_input fin{std::istringstream{input}, seqan3::format_fasta{}};

    data_storage.sequences = std::move(get<seqan3::field::SEQ>(fin)); // we move the buffer directly into our storage
    data_storage.ids       = std::move(get<seqan3::field::ID>(fin));  // we move the buffer directly into our storage
}
