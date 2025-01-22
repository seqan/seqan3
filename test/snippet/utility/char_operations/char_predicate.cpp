// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/utility/char_operations/predicate.hpp>

int main()
{
    //! [is_eof]
    static_assert(seqan3::is_eof(EOF));
    static_assert(!seqan3::is_eof('C'));
    //! [is_eof]

    //! [is_cntrl]
    static_assert(seqan3::is_cntrl('\0'));
    //! [is_cntrl]

    //! [is_print]
    static_assert(seqan3::is_print(' '));
    //! [is_print]

    //! [is_space]
    static_assert(seqan3::is_space('\n'));
    //! [is_space]

    //! [is_blank]
    static_assert(seqan3::is_blank('\t'));
    //! [is_blank]

    //! [is_graph]
    static_assert(seqan3::is_graph('%'));
    //! [is_graph]

    //! [is_punct]
    static_assert(seqan3::is_punct(':'));
    //! [is_punct]

    //! [is_alnum]
    static_assert(seqan3::is_alnum('9'));
    //! [is_alnum]

    //! [is_alpha]
    static_assert(seqan3::is_alpha('z'));
    //! [is_alpha]

    //! [is_upper]
    static_assert(seqan3::is_upper('K'));
    //! [is_upper]

    //! [is_lower]
    static_assert(seqan3::is_lower('a'));
    //! [is_lower]

    //! [is_digit]
    static_assert(seqan3::is_digit('1'));
    //! [is_digit]

    //! [is_xdigit]
    static_assert(seqan3::is_xdigit('e'));
    //! [is_xdigit]
}
