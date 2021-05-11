#include <sstream>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/std/ranges>

auto input = R"(> TEST1
ACGT
> Test2
AGGCTGA
> Test3
GGAGTATAATATATATATATATAT)";

int main()
{
    seqan3::sequence_file_input fin{std::istringstream{}, seqan3::format_fasta{}};

    auto minimum_length5_filter = std::views::filter([] (auto const & rec)
    {
        return std::ranges::size(rec.sequence()) >= 5;
    });

    for (auto & rec : fin | minimum_length5_filter) // only record with sequence length >= 5 will "appear"
    {
        seqan3::debug_stream << "IDs of seq_length >= 5: " << rec.id() << '\n';
        // ...
    }
}
