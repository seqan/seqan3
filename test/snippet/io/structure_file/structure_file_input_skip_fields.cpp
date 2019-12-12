#include <sstream>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/structure_file/input.hpp>
#include <seqan3/range/views/get.hpp>
#include <seqan3/range/views/to_char.hpp>

auto input = R"(> S.cerevisiae_tRNA-PHE M10740/1-73
GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA
(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)
> example
UUGGAGUACACAACCUGUACACUCUUUC
..(((((..(((...)))..)))))... (-3.71))";

int main()
{
    using seqan3::get;

    seqan3::structure_file_input fin{std::istringstream{input},
                                     seqan3::format_vienna{},
                                     seqan3::fields<seqan3::field::id, seqan3::field::structured_seq>{}};

    // note that the order is now different, "id" comes first, because it was specified first
    for (auto & [id, struc_seq] : fin)
    {
        seqan3::debug_stream << "ID: "        << id                                                         << '\n';
        // sequence and structure are part of the same vector, of type std::vector<structured_rna<rna5, wuss51>>
        // sequence and structure strings are extracted and converted to char on-the-fly
        seqan3::debug_stream << "SEQ: "       << (struc_seq | seqan3::views::get<0> | seqan3::views::to_char) << '\n';
        seqan3::debug_stream << "STRUCTURE: " << (struc_seq | seqan3::views::get<1> | seqan3::views::to_char) << '\n';
    }
}
