// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <ranges>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/detail/strong_type.hpp>

struct error : seqan3::detail::strong_type<unsigned, error>
{
    using seqan3::detail::strong_type<unsigned, error>::strong_type;
};

struct window_size : seqan3::detail::strong_type<unsigned, window_size>
{
    using seqan3::detail::strong_type<unsigned, window_size>::strong_type;
};

namespace detail
{

template <std::ranges::forward_range fwd_rng_type>
bool do_find(fwd_rng_type const &, uint8_t const, uint8_t const)
{
    return true;
}

} // namespace detail

template <std::ranges::forward_range fwd_rng_type>
bool search(fwd_rng_type const & rng, window_size const window_size, error const error)
{
    return detail::do_find(rng, window_size.get(), error.get());
}

int main()
{
    using namespace seqan3::literals;

    std::vector<seqan3::dna4> range = "ACGTT"_dna4;
    search(range, window_size{4u}, error{2u});

    return 0;
}
