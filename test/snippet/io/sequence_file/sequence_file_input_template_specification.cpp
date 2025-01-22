// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <sstream>

#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/utility/type_list/type_list.hpp>

// ... input had amino acid sequences
auto input = R"(>TEST1
FQTWE
>Test2
KYRTW
>Test3
EEYQTWEEFARAAEKLYLTDPMKV)";

int main()
{ // Use amino acid traits below
    using sequence_file_input_type =
        seqan3::sequence_file_input<seqan3::sequence_file_input_default_traits_aa,
                                    seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::qual>,
                                    seqan3::type_list<seqan3::format_fasta>>;
    sequence_file_input_type fin{std::istringstream{input}, seqan3::format_fasta{}};
}
