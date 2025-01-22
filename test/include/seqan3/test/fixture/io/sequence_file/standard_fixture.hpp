// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <optional>
#include <tuple>

#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/io/sequence_file/record.hpp>

namespace seqan3::test::fixture::io::sequence_file
{

struct standard_fixture
{
    using types = seqan3::type_list<std::string,                   // seqan3::field::id,
                                    seqan3::dna5_vector,           // seqan3::field::seq,
                                    std::vector<seqan3::phred42>>; // seqan3::field::qual

    using types_as_ids = seqan3::fields<seqan3::field::id, seqan3::field::seq, seqan3::field::qual>;

    using record_type = seqan3::sequence_record<types, types_as_ids>;

    record_type record1{/*.id =*/"ID1",
                        /*.sequence =*/"ACGTTTTTTTTTTTTTTT"_dna5,
                        /*.base_qualities =*/"!##$%&'()*+,-./++-"_phred42};

    record_type record2{
        /*.id =*/"ID2",
        /*.sequence =*/"ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"_dna5,
        /*.base_qualities =*/
        "!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE"_phred42};

    record_type record3{/*.id =*/"ID3 lala",
                        /*.sequence =*/"ACGTTTA"_dna5,
                        /*.base_qualities =*/"!!!!!!!"_phred42};

    std::vector<record_type> records{record1, record2, record3};
};

} // namespace seqan3::test::fixture::io::sequence_file
