#include <sstream>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/structure_file/input.hpp>

auto input = R"(> S.cerevisiae_tRNA-PHE M10740/1-73
GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA
(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)
> example
UUGGAGUACACAACCUGUACACUCUUUC
..(((((..(((...)))..)))))... (-3.71))";

int main()
{
    seqan3::structure_file_input fin{std::istringstream{input}, seqan3::format_vienna{}};
    auto it = seqan3::begin(fin);

    // the following are equivalent:
    auto & rec0 = *it;
    auto & rec1 = fin.front();

    // both become invalid after incrementing "it"!
}
