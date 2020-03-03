#include <sstream>

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

    for (auto & rec : fin)
    {
        seqan3::debug_stream << "ID:  " << seqan3::get<seqan3::field::id>(rec) << '\n';
        seqan3::debug_stream << "SEQ: " << seqan3::get<seqan3::field::seq>(rec) << '\n';
        // a quality field also exists, but is not printed, because we know it's empty for FastA files.
    }
}
