#include <sstream>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/structure_file/input.hpp>
#include <seqan3/range/views/to_char.hpp>
#include <seqan3/std/ranges>

auto input = R"(> S.cerevisiae_tRNA-PHE M10740/1-73
GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA
(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)
> example
UUGGAGUACACAACCUGUACACUCUUUC
..(((((..(((...)))..)))))... (-3.71))";

int main()
{
    using seqan3::get;

    seqan3::structure_file_input fin{std::istringstream{input}, seqan3::format_vienna{}};

    auto minimum_length5_filter = std::views::filter([] (auto const & rec)
    {
        return std::ranges::size(get<seqan3::field::seq>(rec)) >= 5;
    });

    for (auto & rec : fin | minimum_length5_filter) // only record with sequence length >= 5 will "appear"
        seqan3::debug_stream << (get<0>(rec) | seqan3::views::to_char) << '\n';
}
