// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <variant> // for std::visit

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>
#include <seqan3/utility/container/concept.hpp> // for the seqan3::container

// a lambda helper function that prints every type in the std::variant<...allowed SAM tag types...>
auto print_fn = [](auto && arg)
{
    using T = std::remove_cvref_t<decltype(arg)>; // the type T of arg.

    if constexpr (!seqan3::container<T>) // If T is not a container,
    {
        seqan3::debug_stream << arg << '\n'; // just print arg directly.
    }
    else // If T is a container,
    {
        for (auto const & arg_v : arg) // print every value in arg.
            seqan3::debug_stream << arg_v << ",";
        seqan3::debug_stream << '\n';
    }
};

int main()
{
    using namespace seqan3::literals;

    seqan3::sam_tag_dictionary dict{}; // initialise empty dictionary

    // ! there is no get function for unknown tags !
    // dict.get<"XZ"_tag>() = 3;
    // but you can use the operator[]
    dict["XZ"_tag] = 3; // set unknown SAM tag 'XZ' to 3 (type int32_t)

    // ! there is no get function for unknown tags !
    // auto nm = dict.get<"XZ"_tag>();
    // but you can use the operator[] again
    auto xz = dict["XZ"_tag]; // get SAM tag 'XZ' (type std::variant<...allowed SAM tag types...>)

    // ! you cannot print a std::variant directly !
    // seqan3::debug_stream << nm << '\n';
    // but you can use visit:
    std::visit(print_fn, xz); // prints 3
}
