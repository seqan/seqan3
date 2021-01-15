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
    seqan3::structure_file_input fin{std::istringstream{input}, seqan3::format_vienna{}};

#if !SEQAN3_WORKAROUND_GCC_93983
    auto minimum_length5_filter = std::views::filter([] (auto const & rec)
    {
        return std::ranges::size(rec.sequence()) >= 5;
    });
#endif // !SEQAN3_WORKAROUND_GCC_93983

#if SEQAN3_WORKAROUND_GCC_93983
    for (auto & rec : fin /*| minimum_length5_filter*/) // only record with sequence length >= 5 will "appear"
#else // ^^^ workaround / no workaround vvv
    for (auto & rec : fin | minimum_length5_filter) // only record with sequence length >= 5 will "appear"
#endif // SEQAN3_WORKAROUND_GCC_93983
        seqan3::debug_stream << (rec.sequence() | seqan3::views::to_char) << '\n';
}
