// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/mask/mask.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    seqan3::mask my_mask = seqan3::mask::masked;
    seqan3::mask another_mask{};

    my_mask.assign_rank(false);  // will assign my_mask the value mask::unmasked
    another_mask.assign_rank(0); // will also assign another_mask the value mask::unmasked

    if (my_mask.to_rank() == another_mask.to_rank())
        seqan3::debug_stream << "Both are UNMASKED!\n";
}
