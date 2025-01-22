// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/io/detail/record.hpp>
#include <seqan3/utility/type_list/type_list.hpp>

int main()
{
    using types = seqan3::type_list<std::string, seqan3::dna4_vector, std::vector<seqan3::phred42>>;
    using types_as_ids = seqan3::fields<seqan3::field::id, seqan3::field::seq, seqan3::field::qual>;
    using selected_ids = seqan3::fields<seqan3::field::qual, seqan3::field::id>;

    using selected_types = seqan3::detail::select_types_with_ids_t<types, types_as_ids, selected_ids>;
    // resolves to type_list<std::vector<phred42>, std::string>
    static_assert(std::same_as<selected_types, seqan3::type_list<std::vector<seqan3::phred42>, std::string>>);
}
