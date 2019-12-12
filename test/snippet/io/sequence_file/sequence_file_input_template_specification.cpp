#include <sstream>

#include <seqan3/core/type_list/type_list.hpp>
#include <seqan3/io/sequence_file/input.hpp>

// ... input had amino acid sequences
auto input = R"(> TEST1
FQTWE
> Test2
KYRTW
> Test3
EEYQTWEEFARAAEKLYLTDPMKV)";

int main()
{                                                                      // Use amino acid traits below
    using sequence_file_input_type = seqan3::sequence_file_input<seqan3::sequence_file_input_default_traits_aa,
                                                                 seqan3::fields<seqan3::field::seq,
                                                                                seqan3::field::id,
                                                                                seqan3::field::qual>,
                                                                 seqan3::type_list<seqan3::format_fasta>>;
    sequence_file_input_type fin{std::istringstream{input}, seqan3::format_fasta{}};
}
