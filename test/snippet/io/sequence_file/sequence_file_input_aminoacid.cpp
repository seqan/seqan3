#include <sstream>

#include <seqan3/io/sequence_file/input.hpp>

auto input = R"(> TEST1
ACGT
> Test2
AGGCTGA
> Test3
GGAGTATAATATATATATATATAT)";

int main()
{
    seqan3::sequence_file_input<seqan3::sequence_file_input_default_traits_aa> fin{std::istringstream{input},
                                                                                   seqan3::format_fasta{}};
}
