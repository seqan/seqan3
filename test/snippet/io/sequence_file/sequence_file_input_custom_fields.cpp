#include <sstream>
#include <vector>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/range/views/get.hpp>

auto input = R"(> TEST1
ACGT
> Test2
AGGCTGA
> Test3
GGAGTATAATATATATATATATAT)";

int main()
{
    using seqan3::get;

    seqan3::sequence_file_input fin{std::istringstream{input},
                                    seqan3::format_fasta{},
                                    seqan3::fields<seqan3::field::id, seqan3::field::seq_qual>{}};

    for (auto & [id, seq_qual] : fin) // the order is now different, "id" comes first, because it was specified first
    {
        seqan3::debug_stream << "ID:  " << id << '\n';
        // sequence and qualities are part of the same vector, of type std::vector<dna5q>
        seqan3::debug_stream << "SEQ: "  << (seq_qual | seqan3::views::get<0>) << '\n'; // sequence string is extracted
        seqan3::debug_stream << "QUAL: " << (seq_qual | seqan3::views::get<1>) << '\n'; // quality string is extracted
    }
}
