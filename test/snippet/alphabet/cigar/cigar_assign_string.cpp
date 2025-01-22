// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    std::string cigar_str{"4S134M"}; // input

    seqan3::cigar letter1{};
    seqan3::cigar letter2{};

    // Assign from string
    // convenient but creates an unnecessary string copy "4S"
    letter1.assign_string(cigar_str.substr(0, 2));
    letter2.assign_string(cigar_str.substr(2, 4));
    seqan3::debug_stream << letter1 << '\n'; // prints 4S
    seqan3::debug_stream << letter2 << '\n'; // prints 134M

    // Assign from std::string_view (No extra string copies)
    // Version 1
    letter1.assign_string(std::string_view{cigar_str}.substr(0, 2));
    letter2.assign_string(std::string_view{cigar_str}.substr(2, 4));
    seqan3::debug_stream << letter1 << '\n'; // prints 4S
    seqan3::debug_stream << letter2 << '\n'; // prints 134M
    // No extra string copiesersion 2
    letter1.assign_string(/*std::string_view*/ {cigar_str.data(), 2});
    letter2.assign_string(/*std::string_view*/ {cigar_str.data() + 2, 4});
    seqan3::debug_stream << letter1 << '\n'; // prints 4S
    seqan3::debug_stream << letter2 << '\n'; // prints 134M

    // Assign from char array
    letter2.assign_string("40S");
    seqan3::debug_stream << letter2 << '\n'; // prints 40S

    // Assign from seqan3::small_string
    letter2.assign_string(letter1.to_string());
    seqan3::debug_stream << letter2 << '\n'; // prints 4S
}
