// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

//! [exercise]
#include <array>  // std::array
#include <string> // std::string
#include <vector> // std::vector

#include <seqan3/alphabet/all.hpp>
#include <seqan3/alphabet/views/all.hpp> // optional: use views to convert the input string to a dna5 sequence
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/views/all.hpp> // optional: use views to convert the input string to a dna5 sequence

int main(int argc, char * argv[])
{
    using namespace seqan3::literals;

    std::string input{};
    seqan3::argument_parser parser("GC-Content", argc, argv);
    parser.add_positional_option(input, "Specify an input sequence.");
    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext) // the input is invalid
    {
        seqan3::debug_stream << "[PARSER ERROR] " << ext.what() << '\n';
        return 0;
    }

    // Convert the input to a dna5 sequence
    std::vector<seqan3::dna5> sequence{};
    for (char c : input)
        sequence.push_back(seqan3::assign_char_to(c, seqan3::dna5{}));

    // Optional: use views for the conversion. Views will be introduced in the next chapter.
    // std::vector<seqan3::dna5> sequence = input
    //                                    | seqan3::views::char_to<seqan3::dna5>
    //                                    | seqan3::ranges::to<std::vector>();

    // Initialise an array with count values for dna5 symbols.
    std::array<size_t, seqan3::dna5::alphabet_size> count{}; // default initialised with zeroes

    // Increase the symbol count according to the sequence.
    for (seqan3::dna5 symbol : sequence)
        ++count[symbol.to_rank()];

    // Calculate the GC content: (#G + #C) / (#A + #T + #G + #C).
    size_t gc = count['C'_dna5.to_rank()] + count['G'_dna5.to_rank()];
    size_t atgc = input.size() - count['N'_dna5.to_rank()];
    float gc_content = 1.0f * gc / atgc;

    seqan3::debug_stream << "The GC content of " << sequence << " is " << 100 * gc_content << "%.\n";

    return 0;
}
//! [exercise]

void alternatively()
{
    std::string input{};
    // if something changes in here, please update above:
    std::vector<seqan3::dna5> sequence =
        input | seqan3::views::char_to<seqan3::dna5> | seqan3::ranges::to<std::vector>();
}
