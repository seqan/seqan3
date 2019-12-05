#include <sstream>

#include <seqan3/core/type_list/type_list.hpp>
#include <seqan3/io/structure_file/input.hpp>

auto input = R"(> S.cerevisiae_tRNA-PHE M10740/1-73
ACEWACEW
HGEBHHHH
> example
ACEWACEWACEWACEW
HGEBHHHHHGEBHHHH)";

int main()
{
    // ... input had amino acid sequences
    seqan3::structure_file_input<seqan3::structure_file_input_default_traits_aa,
                                 seqan3::fields<seqan3::field::SEQ, seqan3::field::ID, seqan3::field::STRUCTURE>,
                                 seqan3::type_list<seqan3::format_vienna>> fin{std::istringstream{input},
                                                                               seqan3::format_vienna{}};
}
