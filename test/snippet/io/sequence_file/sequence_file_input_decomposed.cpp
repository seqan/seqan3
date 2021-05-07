#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>

auto input = R"(> TEST1
ACGT
> Test2
AGGCTGA
> Test3
GGAGTATAATATATATATATATAT)";

int main()
{
    seqan3::sequence_file_input fin{std::istringstream{input}, seqan3::format_fasta{}};

    for (auto & [sequence, id, quality] : fin)
    {
        seqan3::debug_stream << "ID: "        << id       << '\n';
        seqan3::debug_stream << "SEQ: "       << sequence << '\n';
        seqan3::debug_stream << "EMPTY QUAL." << quality  << '\n'; // quality is empty for FastA files
    }
}
