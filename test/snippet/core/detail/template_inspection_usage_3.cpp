// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <vector>

#include <seqan3/core/detail/template_inspection.hpp>

int main()
{
    using my_type = std::vector<int>;

    if constexpr (seqan3::detail::template_specialisation_of<my_type, std::vector>) // Note: std::vector has no <> !
    {
        // ...
    }
}
