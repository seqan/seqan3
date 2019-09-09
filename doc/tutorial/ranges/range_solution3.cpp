#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/views/all.hpp>            // include all of SeqAn's views
#include <seqan3/std/ranges>                    // include all of the standard library's views

int main(int argc, char** argv)
{
    if (argc == 1)
        return 0;

    std::string s{argv[1]};

    auto s_as_dna = s | seqan3::views::char_to<seqan3::dna5>;
    // Bonus:
    //auto s_as_dna = s | std::views::transform([] (char const c)
    //{
    //    return seqan3::assign_char_strictly_to(c, seqan3::dna5{});
    //});

    seqan3::debug_stream << "Original: " << s_as_dna << '\n';
    seqan3::debug_stream << "RevComp:  " << (s_as_dna | std::views::reverse | seqan3::views::complement) << '\n';
    seqan3::debug_stream << "Frames:   " << (s_as_dna | seqan3::views::translate) << '\n';
}
