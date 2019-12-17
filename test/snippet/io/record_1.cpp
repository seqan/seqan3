#include <sstream>

#include <seqan3/io/record.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>

auto input = R"(> TEST1
ACGT
> Test2
AGGCTGA
> Test3
GGAGTATAATATATATATATATAT)";

int main()
{
    using seqan3::get;

    // specify custom field combination/order to file:
    seqan3::sequence_file_input fin{std::istringstream{input},
                                    seqan3::format_fasta{},
                                    seqan3::fields<seqan3::field::id, seqan3::field::seq>{}};

    auto record = fin.front(); // get current record, in this case the first

    // record is tuple-like type, but allows access via field identifiers:
    auto & id = get<seqan3::field::id>(record);
    auto & seq = get<seqan3::field::seq>(record);
}
