#include <sstream>
#include <utility>
#include <vector>

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

    using record_type = typename decltype(fin)::record_type;
    std::vector<record_type> records;

    for (auto & rec : fin)
        records.push_back(std::move(rec));
}
