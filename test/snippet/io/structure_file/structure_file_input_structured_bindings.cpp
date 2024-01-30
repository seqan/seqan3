// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <sstream>

#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/structure_file/input.hpp>
#include <seqan3/utility/type_list/type_list.hpp>

auto input = R"(> S.cerevisiae_tRNA-PHE M10740/1-73
ACEWACEW
HGEBHHHH
> example
ACEWACEWACEWACEW
HGEBHHHHHGEBHHHH)";

int main()
{
    using seqan3::get;

    using structure_file_input_t =
        seqan3::structure_file_input<seqan3::structure_file_input_default_traits_aa,
                                     seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::structure>,
                                     seqan3::type_list<seqan3::format_vienna>>;

    structure_file_input_t fin{std::istringstream{input}, seqan3::format_vienna{}};

    for (auto & [seq, id, structure] : fin)
    {
        seqan3::debug_stream << "ID: " << id << '\n';
        // sequence and structure are converted to char on-the-fly
        seqan3::debug_stream << "SEQ: " << (seq | seqan3::views::to_char) << '\n';
        seqan3::debug_stream << "STRUCTURE: " << (structure | seqan3::views::to_char) << '\n';
    }
}
