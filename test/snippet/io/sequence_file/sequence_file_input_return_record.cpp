#include <sstream>

#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/std/ranges>

auto input = R"(> TEST1
ACGT
> Test2
AGGCTGA
> Test3
GGAGTATAATATATATATATATAT)";

int main()
{
    seqan3::sequence_file_input fin{std::istringstream{input}, seqan3::format_fasta{}};
    auto it = std::ranges::begin(fin);

    // the following are equivalent:
    auto & rec0 = *it;
    auto & rec1 = fin.front();

    // Note: rec0 and rec1 are references and become invalid after incrementing "it"!
}
