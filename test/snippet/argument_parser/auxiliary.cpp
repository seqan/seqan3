// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/argument_parser/all.hpp>

int main(int argc, char const ** argv)
{
    seqan3::argument_parser myparser{"Test", argc, argv};
    std::string myvar{"Example"};
    myparser.add_option(myvar, 's', "special-op", "You know what you doin'?", seqan3::option_spec::advanced);
}
