// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <string>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/type_list/detail/type_list_algorithm.hpp>
#include <seqan3/utility/type_list/type_list.hpp>

namespace incomplete
{

struct type;

} // namespace incomplete

int main()
{
    // With c++20 you could also write it like this
    // auto fn = []<typename type>(std::type_identity<type>)
    // {
    // ...
    // };
    auto fn = [](auto id)
    {
        // id is of type std::type_identity<type>
        using id_t = decltype(id);
        using type = typename id_t::type;

        static_assert(std::is_same_v<id_t, std::type_identity<type>>, "id is of type std::type_identity<type>");

        if constexpr (std::is_same_v<type, bool>)
            seqan3::debug_stream << "bool";
        else if constexpr (std::is_same_v<type, int>)
            seqan3::debug_stream << "int";
        else if constexpr (std::is_same_v<type, float>)
            seqan3::debug_stream << "float";
        else if constexpr (std::is_same_v<type, incomplete::type>)
            seqan3::debug_stream << "incomplete::type";

        seqan3::debug_stream << ", ";
    };

    // prints each type name, i.e. "int, float, bool, incomplete::type, \n"
    using types = seqan3::type_list<int, float, bool, incomplete::type>;
    seqan3::detail::for_each<types>(fn);
    seqan3::debug_stream << "\n";

    // is the same as explicitly writing
    fn(std::type_identity<int>{});
    fn(std::type_identity<float>{});
    fn(std::type_identity<bool>{});
    fn(std::type_identity<incomplete::type>{});
    seqan3::debug_stream << "\n";
    return 0;
}
