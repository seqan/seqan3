// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <concepts>
#include <ranges>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/pretty_printing.hpp>
#include <seqan3/utility/views/convert.hpp>
#include <seqan3/utility/views/zip.hpp>

using seqan3::operator""_dna5;
using seqan3::operator""_phred42;

struct sequence_file_data : public ::testing::Test
{
    std::vector<std::string> ids{
        {"ID1"},
        {"ID2"},
        {"ID3 lala"},
    };

    std::vector<seqan3::dna5_vector> seqs{
        {"ACGTTTTTTTTTTTTTTT"_dna5},
        {"ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"_dna5},
        {"ACGTTTA"_dna5},
    };

    std::vector<std::vector<seqan3::phred42>> quals{
        {"!##$%&'()*+,-./++-"_phred42},
        {"!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE"_phred42},
        {"!!!!!!!"_phred42},
    };

    std::ostringstream ostream{};
};

template <typename format_t>
struct sequence_file_read : public sequence_file_data
{};

TYPED_TEST_SUITE_P(sequence_file_read);

// ----------------------------------------------------------------------------
// general
// ----------------------------------------------------------------------------

TYPED_TEST_P(sequence_file_read, concept_check)
{
    EXPECT_TRUE((seqan3::sequence_file_input_format<TypeParam>));
}

// ----------------------------------------------------------------------------
// sequence_file_read
// ----------------------------------------------------------------------------

TYPED_TEST_P(sequence_file_read, standard)
{
    std::stringstream istream{this->standard_input};
    seqan3::sequence_file_input fin{istream, TypeParam{}};

    auto it = fin.begin();
    for (unsigned i = 0; i < 3; ++i, ++it)
    {
        EXPECT_RANGE_EQ((*it).sequence(), this->seqs[i]);
        EXPECT_EQ((*it).id(), this->ids[i]);
        if constexpr (std::same_as<TypeParam, seqan3::format_fastq> || std::same_as<TypeParam, seqan3::format_sam>)
        {
            EXPECT_RANGE_EQ((*it).base_qualities(), this->quals[i]);
        }
    }
}

TYPED_TEST_P(sequence_file_read, only_seq)
{
    std::stringstream istream{this->standard_input};
    seqan3::sequence_file_input fin{istream, TypeParam{}, seqan3::fields<seqan3::field::seq>{}};

    auto it = fin.begin();
    for (unsigned i = 0; i < 3; ++i, ++it)
        EXPECT_EQ(std::get<0>(*it), this->seqs[i]);
}

TYPED_TEST_P(sequence_file_read, only_id)
{
    std::stringstream istream{this->standard_input};
    seqan3::sequence_file_input fin{istream, TypeParam{}, seqan3::fields<seqan3::field::id>{}};

    auto it = fin.begin();
    for (unsigned i = 0; i < 3; ++i, ++it)
        EXPECT_EQ(std::get<0>(*it), this->ids[i]);
}

TYPED_TEST_P(sequence_file_read, options_truncate_ids)
{
    std::stringstream istream{this->standard_input};
    seqan3::sequence_file_input fin{istream, TypeParam{}, seqan3::fields<seqan3::field::id>{}};
    fin.options.truncate_ids = true;
    this->ids[2] = "ID3"; // "lala" is stripped

    auto it = fin.begin();
    for (unsigned i = 0; i < 3; ++i, ++it)
        EXPECT_EQ(std::get<0>(*it), this->ids[i]);
}

TYPED_TEST_P(sequence_file_read, illegal_alphabet_character)
{
    std::stringstream istream{this->illegal_alphabet_character_input};
    seqan3::sequence_file_input fin{istream, TypeParam{}};
    EXPECT_THROW(fin.begin(), seqan3::parse_error);
}

TYPED_TEST_P(sequence_file_read, no_or_ill_formatted_id)
{
    std::stringstream istream{this->no_or_ill_formatted_id_input};
    seqan3::sequence_file_input fin{istream, TypeParam{}};
    EXPECT_THROW(fin.begin(), seqan3::parse_error);
}

// ----------------------------------------------------------------------------
// sequence_file_write
// ----------------------------------------------------------------------------

template <typename format_type>
struct sequence_file_write : public sequence_file_read<format_type>
{};

TYPED_TEST_SUITE_P(sequence_file_write);

TYPED_TEST_P(sequence_file_write, concept_check)
{
    EXPECT_TRUE((seqan3::sequence_file_output_format<TypeParam>));
}

TYPED_TEST_P(sequence_file_write, standard)
{
    seqan3::sequence_file_output fout{this->ostream, TypeParam{}};

    for (unsigned i = 0; i < 3; ++i)
        EXPECT_NO_THROW((fout.emplace_back(this->seqs[i], this->ids[i], this->quals[i])));

    this->ostream.flush();
    EXPECT_EQ(this->ostream.str(), this->standard_output);
}

TYPED_TEST_P(sequence_file_write, arg_handling_id_missing)
{
    if constexpr (!std::same_as<TypeParam, seqan3::format_sam>)
    {
        seqan3::sequence_file_output fout{this->ostream, TypeParam{}, seqan3::fields<seqan3::field::seq>{}};
        EXPECT_THROW((fout.emplace_back(this->seqs[0])), std::logic_error);
    }
}

TYPED_TEST_P(sequence_file_write, arg_handling_id_empty)
{
    if constexpr (!std::same_as<TypeParam, seqan3::format_sam>)
    {
        seqan3::sequence_file_output fout{this->ostream,
                                          TypeParam{},
                                          seqan3::fields<seqan3::field::seq, seqan3::field::id>{}};
        EXPECT_THROW((fout.emplace_back(this->seqs[0], std::string_view{""}, std::ignore)), std::runtime_error);
    }
}

TYPED_TEST_P(sequence_file_write, arg_handling_seq_missing)
{
    if constexpr (!std::same_as<TypeParam, seqan3::format_sam>)
    {
        seqan3::sequence_file_output fout{this->ostream, TypeParam{}, seqan3::fields<seqan3::field::id>{}};
        EXPECT_THROW((fout.emplace_back(this->ids[0])), std::logic_error);
    }
}

TYPED_TEST_P(sequence_file_write, arg_handling_seq_empty)
{
    if constexpr (!std::same_as<TypeParam, seqan3::format_sam>)
    {
        seqan3::sequence_file_output fout{this->ostream,
                                          TypeParam{},
                                          seqan3::fields<seqan3::field::seq, seqan3::field::id>{}};
        EXPECT_THROW((fout.emplace_back(std::string_view{""}, this->ids[0], std::ignore)), std::runtime_error);
    }
}

REGISTER_TYPED_TEST_SUITE_P(sequence_file_read,
                            concept_check,
                            standard,
                            only_seq,
                            only_id,
                            illegal_alphabet_character,
                            options_truncate_ids,
                            no_or_ill_formatted_id);

REGISTER_TYPED_TEST_SUITE_P(sequence_file_write,
                            concept_check,
                            standard,
                            arg_handling_id_missing,
                            arg_handling_id_empty,
                            arg_handling_seq_missing,
                            arg_handling_seq_empty);
