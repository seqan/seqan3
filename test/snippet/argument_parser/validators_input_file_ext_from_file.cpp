// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/argument_parser/validators.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>

int main()
{
    // Default constructed validator has an empty extension list.
    seqan3::input_file_validator validator1{};
    seqan3::debug_stream << validator1.get_help_page_message() << '\n';

    // Specify your own extensions for the input file.
    seqan3::input_file_validator validator2{std::vector{std::string{"exe"}, std::string{"fasta"}}};
    seqan3::debug_stream << validator2.get_help_page_message() << '\n';

    // Give the seqan3 file type as a template argument to get all valid extensions for this file.
    seqan3::input_file_validator<seqan3::sequence_file_input<>> validator3{};
    seqan3::debug_stream << validator3.get_help_page_message() << '\n';

    return 0;
}
