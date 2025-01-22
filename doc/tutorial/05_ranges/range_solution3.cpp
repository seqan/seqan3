// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <ranges> // include all of the standard library's views

#include <seqan3/alphabet/views/all.hpp>  // include all of SeqAn's views
#include <seqan3/argument_parser/all.hpp> // optional: include the argument_parser
#include <seqan3/core/debug_stream.hpp>

int main(int argc, char ** argv)
{
    // We use the seqan3::argument_parser which was introduced in the second chapter
    // of the tutorial: "Parsing command line arguments with SeqAn".
    seqan3::argument_parser myparser{"Assignment-3", argc, argv}; // initialize
    std::string s{};

    myparser.add_positional_option(s, "Please specify the DNA string.");

    try
    {
        myparser.parse();
    }
    catch (seqan3::argument_parser_error const & ext) // the user did something wrong
    {
        std::cerr << "[PARSER ERROR]" << ext.what() << '\n'; // you can customize your error message
        return 0;
    }

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
