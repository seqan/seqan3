// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>

#include "../alignment/aligned_sequence_test_template.hpp"

template <typename container_type>
    requires seqan3::aligned_sequence<container_type>
class aligned_sequence<container_type> : public ::testing::Test
{
public:
    // Initializer function is needed for the typed test because the gapped_decorator
    // will be initialized differently than the naive std::vector<seqan3::gapped<dna>>.
    void initialise_typed_test_container(container_type & container, seqan3::dna4_vector const & target)
    {
        container.clear();
        for (auto & val : target)
        {
            container.push_back(seqan3::assign_char_to(seqan3::to_char(val), typename container_type::value_type{}));
        }
    }
};

using test_types = ::testing::Types<std::vector<seqan3::gapped<seqan3::dna4>>,
                                    std::vector<seqan3::gapped<seqan3::qualified<seqan3::dna4, seqan3::phred42>>>>;

INSTANTIATE_TYPED_TEST_SUITE_P(container_of_gapped_alphabets, aligned_sequence, test_types, );
