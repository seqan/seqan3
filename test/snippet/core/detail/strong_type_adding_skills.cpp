// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/core/detail/strong_type.hpp>

using seqan3::operator|;

struct error :
    seqan3::detail::strong_type<uint8_t,
                                error,
                                seqan3::detail::strong_type_skill::decrement
                                    | seqan3::detail::strong_type_skill::increment>
{
    using seqan3::detail::strong_type<uint8_t,
                                      error,
                                      seqan3::detail::strong_type_skill::decrement
                                          | seqan3::detail::strong_type_skill::increment>::strong_type;
};

int main()
{
    error e{4u};
    --e;
    ++e;
}
