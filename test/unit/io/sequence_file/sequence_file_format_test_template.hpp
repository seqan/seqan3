// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/range/views/convert.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/std/ranges>
#include <seqan3/std/concepts>
#include <seqan3/test/pretty_printing.hpp>

using seqan3::operator""_dna5;
using seqan3::operator""_phred42;

struct sequence_file_data : public ::testing::Test
{
    std::vector<std::string> ids
    {
        { "ID1" },
        { "ID2" },
        { "ID3 lala" },
    };

    std::vector<seqan3::dna5_vector> seqs
    {
        { "ACGTTTTTTTTTTTTTTT"_dna5 },
        { "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"_dna5 },
        { "ACGTTTA"_dna5 },
    };

    std::vector<std::vector<seqan3::phred42>> quals
    {
        { "!##$%&'()*+,-./++-"_phred42 },
        { "!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE"_phred42 },
        { "!!!!!!!"_phred42 },
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
        EXPECT_TRUE((std::ranges::equal(seqan3::get<seqan3::field::seq>(*it), this->seqs[i])));
        EXPECT_EQ(seqan3::get<seqan3::field::id>(*it), this->ids[i]);
        if constexpr (std::same_as<TypeParam, seqan3::format_fastq> || std::same_as<TypeParam, seqan3::format_sam>)
        {
            EXPECT_TRUE((std::ranges::equal(seqan3::get<seqan3::field::qual>(*it), this->quals[i])));
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

TYPED_TEST_P(sequence_file_read, seq_qual)
{
    std::stringstream istream{this->standard_input};
    seqan3::sequence_file_input fin{istream, TypeParam{}, seqan3::fields<seqan3::field::id, seqan3::field::seq_qual>{}};

    auto it = fin.begin();
    for (unsigned i = 0; i < 3; ++i, ++it)
    {
        EXPECT_TRUE((std::ranges::equal(seqan3::get<seqan3::field::id>(*it), this->ids[i])));
        EXPECT_TRUE((std::ranges::equal(seqan3::get<seqan3::field::seq_qual>(*it)
                                        | seqan3::views::convert<seqan3::dna5>, this->seqs[i])));
    }
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

TYPED_TEST_P(sequence_file_write, seq_qual)
{
    auto convert_to_qualified = std::views::transform([] (auto const in)
    {
        return seqan3::qualified<seqan3::dna5, seqan3::phred42>{std::get<0>(in), std::get<1>(in)};
    });

    seqan3::sequence_file_output fout{this->ostream, TypeParam{}, seqan3::fields<seqan3::field::id,
                                                                                 seqan3::field::seq_qual>{}};

    for (unsigned i = 0; i < 3; ++i)
    {
        EXPECT_NO_THROW((fout.emplace_back(this->ids[i],
                                           seqan3::views::zip(this->seqs[i], this->quals[i]) | convert_to_qualified)));
    }

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
        seqan3::sequence_file_output fout{this->ostream, TypeParam{}, seqan3::fields<seqan3::field::seq,
                                                                                     seqan3::field::id>{}};
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
        seqan3::sequence_file_output fout{this->ostream, TypeParam{}, seqan3::fields<seqan3::field::seq,
                                                                                     seqan3::field::id>{}};
        EXPECT_THROW((fout.emplace_back(std::string_view{""}, this->ids[0], std::ignore)), std::runtime_error);
    }
}

REGISTER_TYPED_TEST_SUITE_P(sequence_file_read,
                            concept_check,
                            standard,
                            only_seq,
                            only_id,
                            seq_qual,
                            illegal_alphabet_character,
                            options_truncate_ids,
                            no_or_ill_formatted_id);

REGISTER_TYPED_TEST_SUITE_P(sequence_file_write,
                            concept_check,
                            standard,
                            seq_qual,
                            arg_handling_id_missing,
                            arg_handling_id_empty,
                            arg_handling_seq_missing,
                            arg_handling_seq_empty);
