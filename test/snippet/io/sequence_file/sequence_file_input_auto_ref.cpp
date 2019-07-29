#include <sstream>
#include <utility>
#include <vector>

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

    using record_type = typename decltype(fin)::record_type;
    std::vector<record_type> records{};

    for (auto & rec : fin)
        records.push_back(std::move(rec));
}
