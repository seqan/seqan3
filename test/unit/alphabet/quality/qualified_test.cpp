// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/quality/aliases.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/phred63.hpp>
#include <seqan3/alphabet/quality/phred94.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>

#include "../alphabet_constexpr_test_template.hpp"
#include "../alphabet_test_template.hpp"
#include "../composite/alphabet_tuple_base_test_template.hpp"
#include "../semi_alphabet_constexpr_test_template.hpp"
#include "../semi_alphabet_test_template.hpp"

template <typename alphabet_type, typename phred_type>
class alphabet_tuple_base_test<seqan3::qualified<alphabet_type, phred_type>> : public ::testing::Test
{
public:
    using T = seqan3::qualified<alphabet_type, phred_type>;

    using other_type = std::conditional_t<
        std::is_same_v<alphabet_type, seqan3::dna4>,
        seqan3::rna4,
        std::conditional_t<std::is_same_v<alphabet_type, seqan3::aa27>,
                           seqan3::aa27,
                           std::conditional_t<std::is_same_v<alphabet_type, seqan3::gapped<seqan3::dna4>>,
                                              seqan3::gapped<seqan3::dna4>,
                                              alphabet_type>>>;

    T instance = T{value_1(), value_2()};
    T zero_instance = T{decltype(value_1()){}, decltype(value_2()){}};
    size_t tup_size{2};

    // structured_rna<alphabet_type, phred_type>
    // -------------------------------------------------------------------------
    alphabet_type value_1()
    {
        return alphabet_type{}.assign_char('G');
    }
    other_type assignable_to_value_1()
    {
        return other_type{}.assign_char('G');
    }
    phred_type value_2()
    {
        return phred_type{}.assign_phred(6);
    }
    phred_type assignable_to_value_2()
    {
        return phred_type{}.assign_phred(6);
    }
    auto values_to_cmp()
    {
        return std::make_tuple(/*low */ alphabet_type{}.assign_char('A'),
                               phred_type{}.assign_phred(1),
                               /*mid */ alphabet_type{}.assign_char('C'),
                               phred_type{}.assign_phred(4),
                               /*high*/ alphabet_type{}.assign_char('T'),
                               phred_type{}.assign_phred(9));
    }
};

using qualified_types = ::testing::Types<seqan3::qualified<seqan3::dna4, seqan3::phred42>,
                                         seqan3::qualified<seqan3::dna4, seqan3::phred63>,
                                         seqan3::qualified<seqan3::dna4, seqan3::phred94>,
                                         seqan3::qualified<seqan3::aa27, seqan3::phred42>,
                                         seqan3::qualified<seqan3::gapped<seqan3::dna4>, seqan3::phred42>,
                                         seqan3::dna4q>;

INSTANTIATE_TYPED_TEST_SUITE_P(qualified, alphabet, qualified_types, );
INSTANTIATE_TYPED_TEST_SUITE_P(qualified, semi_alphabet_test, qualified_types, );
INSTANTIATE_TYPED_TEST_SUITE_P(qualified, alphabet_constexpr, qualified_types, );
INSTANTIATE_TYPED_TEST_SUITE_P(qualified, semi_alphabet_constexpr, qualified_types, );
INSTANTIATE_TYPED_TEST_SUITE_P(qualified, alphabet_tuple_base_test, qualified_types, );
