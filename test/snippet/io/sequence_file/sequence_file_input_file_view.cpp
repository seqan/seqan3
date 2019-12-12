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
    using seqan3::get;

    seqan3::sequence_file_input fin{std::istringstream{}, seqan3::format_fasta{}};

    auto minimum_length5_filter = std::views::filter([] (auto const & rec)
    {
        return std::ranges::size(get<seqan3::field::seq>(rec)) >= 5;
    });

    for (auto & rec : fin | minimum_length5_filter) // only record with sequence length >= 5 will "appear"
    {
        seqan3::debug_stream << "IDs of seq_length >= 5: " << get<seqan3::field::id>(rec) << '\n';
        // ...
    }
}
