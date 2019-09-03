#include <sstream>
#include <utility>

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
    //                         ^ no need to specify the template arguments
}
