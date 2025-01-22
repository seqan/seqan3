// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <algorithm>
#include <iterator>
#include <set>
#include <sstream>
#include <tuple>
#include <vector>

#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/nucleotide/rna15.hpp>
#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/alphabet/structure/structured_rna.hpp>
#include <seqan3/alphabet/structure/wuss.hpp>
#include <seqan3/io/structure_file/format_vienna.hpp>
#include <seqan3/io/structure_file/input.hpp>
#include <seqan3/io/structure_file/output.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/utility/views/convert.hpp>

using seqan3::operator""_rna5;
using seqan3::operator""_wuss51;

// ----------------------------------------------------------------------------
// general
// ----------------------------------------------------------------------------

TEST(general, concepts)
{
    EXPECT_TRUE((seqan3::structure_file_input_format<seqan3::format_vienna>));
    EXPECT_TRUE((seqan3::structure_file_output_format<seqan3::format_vienna>));
}

// ----------------------------------------------------------------------------
// reading
// ----------------------------------------------------------------------------

struct read : public ::testing::Test
{
    std::string input{"> S.cerevisiae_tRNA-PHE M10740/1-73\n"
                      "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
                      "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)\n"
                      "> example 2\n"
                      "UUGGAGUACACAACCUGUACACUCUUUC\n"
                      "..(((((..(((...)))..)))))... (-3.71)\n"};

    std::vector<std::string> expected_id{{"S.cerevisiae_tRNA-PHE M10740/1-73"}, {"example 2"}};

    std::vector<seqan3::rna5_vector> const expected_seq{
        {"GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA"_rna5},
        {"UUGGAGUACACAACCUGUACACUCUUUC"_rna5}};

    std::vector<std::vector<seqan3::wuss<51>>> const expected_structure{
        {"(((((((..((((........)))).((((.........)))).....(((((.......))))))))))))."_wuss51},
        {"..(((((..(((...)))..)))))..."_wuss51}};

    std::vector<double> const expected_energy{-17.5, -3.71};

    std::vector<std::vector<uint8_t>> const expected_interactions{
        {71, 70, 69, 68, 67, 66, 65, 24, 23, 22, 21, 12, 11, 10, 9, 42, 41, 40, 39, 29,
         28, 27, 26, 64, 63, 62, 61, 60, 52, 51, 50, 49, 48, 6,  5, 4,  3,  2,  1,  0},
        {24, 23, 22, 21, 20, 17, 16, 15, 11, 10, 9, 6, 5, 4, 3, 2}};

    seqan3::structure_file_input_options<seqan3::rna15, false> options;

    bool check_seq = true;
    bool check_id = true;
    bool check_structure = true;
    bool check_energy = true;

    void bpp_test(std::vector<std::set<std::pair<double, size_t>>> & bpp, std::vector<uint8_t> const & bpp_comp)
    {
        size_t cnt = 0ul;
        auto interaction_sets = bpp
                              | std::views::filter(
                                    [](auto & set)
                                    {
                                        return set.size() == 1;
                                    });
        for (auto & iset : interaction_sets)
        {
            EXPECT_EQ(iset.size(), 1u);
            EXPECT_EQ(iset.begin()->second, bpp_comp[cnt++]);
            EXPECT_FLOAT_EQ(iset.begin()->first, 1.0f);
        }
        EXPECT_EQ(cnt, bpp_comp.size());
    }

    void do_read_test(std::string const & input)
    {
        std::stringstream istream{input};

        auto field_set = seqan3::fields<seqan3::field::id,
                                        seqan3::field::seq,
                                        seqan3::field::bpp,
                                        seqan3::field::structure,
                                        seqan3::field::energy>{};
        seqan3::structure_file_input fin{istream, seqan3::format_vienna{}, field_set};
        fin.options = options;

        auto it = fin.begin();
        for (size_t idx = 0ul; idx < expected_seq.size(); ++idx, ++it)
        {
            EXPECT_EQ(check_energy, (*it).energy().has_value());
            if (check_seq)
            {
                EXPECT_RANGE_EQ((*it).sequence(), expected_seq[idx]);
            }
            if (check_id)
            {
                EXPECT_RANGE_EQ((*it).id(), expected_id[idx]);
            }
            if (check_structure)
            {
                bpp_test((*it).base_pair_probability_matrix(), expected_interactions[idx]);
            }
            if (check_energy)
            {
                EXPECT_DOUBLE_EQ(*(*it).energy(), expected_energy[idx]);
            }
            if (check_structure)
            {
                EXPECT_RANGE_EQ((*it).sequence_structure(), expected_structure[idx]);
            }
        }
    }
};

TEST_F(read, standard)
{
    do_read_test(input);
}

TEST_F(read, newline_before_eof)
{
    input = "> S.cerevisiae_tRNA-PHE M10740/1-73\n"
            "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
            "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)\n"
            "> example 2\n"
            "UUGGAGUACACAACCUGUACACUCUUUC\n"
            "..(((((..(((...)))..)))))... (-3.71)";
    do_read_test(input);
}

TEST_F(read, whitespace_in_seq)
{
    input = "> S.cerevisiae_tRNA-PHE M10740/1-73\n"
            "GCGGAUUUAGCU CAGUUGGGAGAGCGCCAGACU GAAGAUUUGGAGGUC CUGUGUUCGAUCCACA   GAAUU CGCA\n"
            "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)\n"
            "> example 2\n"
            "UUGGAGUAC   ACAACCUGUACAC UCUUUC \n"
            "..(((((..(((...)))..)))))... (-3.71)\n";
    do_read_test(input);
}

TEST_F(read, no_energies)
{
    check_energy = false;
    input = "> S.cerevisiae_tRNA-PHE M10740/1-73\n"
            "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
            "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))).\n"
            "> example 2\n"
            "UUGGAGUACACAACCUGUACACUCUUUC\n"
            "..(((((..(((...)))..)))))...\n";
    do_read_test(input);
}

TEST_F(read, no_ids)
{
    check_id = false;
    input = "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
            "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)\n"
            "UUGGAGUACACAACCUGUACACUCUUUC\n"
            "..(((((..(((...)))..)))))... (-3.71)\n";
    do_read_test(input);
}

TEST_F(read, spaces_and_carriage_return)
{
    check_id = false;
    input = "> S.cerevisiae_tRNA-PHE M10740/1-73\r\n"
            "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\r\n"
            "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.5)\r\n"
            ">example 2\r\n"
            "UUGGAGUA CACAACCUGUACA  CUCU UUC \r\n"
            "..(((((..(((...)))..)))))...     ( -3.71 )\r\n";
    do_read_test(input);
}

TEST_F(read, options_truncate_ids)
{
    input = "> S.cerevisiae_tRNA-PHE M10740/1-73\n"
            "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
            "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)\n"
            "> example 2\n"
            "UUGGAGUACACAACCUGUACACUCUUUC\n"
            "..(((((..(((...)))..)))))... (-3.71)\n";

    options.truncate_ids = true;
    expected_id = {"S.cerevisiae_tRNA-PHE", "example"};
    do_read_test(input);
}

struct read_fields : public read
{};

TEST_F(read_fields, only_seq)
{
    std::stringstream istream{input};
    seqan3::structure_file_input fin{istream, seqan3::format_vienna{}, seqan3::fields<seqan3::field::seq>{}};
    auto it = fin.begin();
    for (size_t idx = 0ul; idx < expected_seq.size(); ++idx, ++it)
    {
        EXPECT_RANGE_EQ((*it).sequence(), expected_seq[idx]);
    }
}

TEST_F(read_fields, only_id)
{
    std::stringstream istream{input};
    seqan3::structure_file_input fin{istream, seqan3::format_vienna{}, seqan3::fields<seqan3::field::id>{}};
    auto it = fin.begin();
    for (size_t idx = 0ul; idx < expected_seq.size(); ++idx, ++it)
    {
        EXPECT_RANGE_EQ((*it).id(), expected_id[idx]);
    }
}

TEST_F(read_fields, only_structure)
{
    std::stringstream istream{input};
    seqan3::structure_file_input fin{istream, seqan3::format_vienna{}, seqan3::fields<seqan3::field::structure>{}};
    auto it = fin.begin();
    for (size_t idx = 0ul; idx < expected_seq.size(); ++idx, ++it)
    {
        EXPECT_RANGE_EQ((*it).sequence_structure(), expected_structure[idx]);
    }
}

TEST_F(read_fields, only_energy)
{
    std::stringstream istream{input};
    seqan3::structure_file_input fin{istream, seqan3::format_vienna{}, seqan3::fields<seqan3::field::energy>{}};
    auto it = fin.begin();
    for (size_t idx = 0ul; idx < expected_seq.size(); ++idx, ++it)
    {
        EXPECT_TRUE((*it).energy());
        EXPECT_DOUBLE_EQ(*(*it).energy(), expected_energy[idx]);
    }
}

TEST_F(read_fields, only_bpp)
{
    std::stringstream istream{input};
    seqan3::structure_file_input fin{istream, seqan3::format_vienna{}, seqan3::fields<seqan3::field::bpp>{}};
    auto it = fin.begin();
    for (size_t idx = 0ul; idx < expected_seq.size(); ++idx, ++it)
    {
        bpp_test((*it).base_pair_probability_matrix(), expected_interactions[idx]);
    }
}

struct read_fail : public read
{};

TEST_F(read_fail, wrong_id)
{
    input[0] = '#'; // invalid character for ID line
    std::stringstream istream{input};
    seqan3::structure_file_input fin{istream, seqan3::format_vienna{}};
    EXPECT_THROW(fin.begin(), seqan3::parse_error);
}

TEST_F(read_fail, missing_seq)
{
    input = "> S.cerevisiae_tRNA-PHE M10740/1-73\n"
            "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)\n";
    std::stringstream istream{input};
    seqan3::structure_file_input fin{istream, seqan3::format_vienna{}};
    EXPECT_THROW(fin.begin(), seqan3::parse_error);
}

TEST_F(read_fail, missing_structure)
{
    input = "> S.cerevisiae_tRNA-PHE M10740/1-73\n"
            "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
            "> example 2\n"
            "UUGGAGUACACAACCUGUACACUCUUUC\n";
    std::stringstream istream{input};
    seqan3::structure_file_input fin{istream, seqan3::format_vienna{}};
    EXPECT_THROW(fin.begin(), seqan3::parse_error);
}

TEST_F(read_fail, missing_structure_and_id)
{
    input = "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
            "UUGGAGUACACAACCUGUACACUCUUUC\n";
    std::stringstream istream{input};
    seqan3::structure_file_input fin{istream, seqan3::format_vienna{}};
    EXPECT_THROW(fin.begin(), seqan3::parse_error);
}

TEST_F(read_fail, structure_too_long)
{
    input = "> S.cerevisiae_tRNA-PHE M10740/1-73\n"
            "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
            "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))).. (-17.50)\n";
    std::stringstream istream{input};
    seqan3::structure_file_input fin{istream, seqan3::format_vienna{}};
    EXPECT_THROW(fin.begin(), seqan3::parse_error);
}

TEST_F(read_fail, structure_too_short)
{
    input = "> S.cerevisiae_tRNA-PHE M10740/1-73\n"
            "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
            "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))) (-17.50)\n";
    std::stringstream istream{input};
    seqan3::structure_file_input fin{istream, seqan3::format_vienna{}};
    EXPECT_THROW(fin.begin(), seqan3::parse_error);
}

TEST_F(read_fail, structure_too_long_structured_seq)
{
    input = "> S.cerevisiae_tRNA-PHE M10740/1-73\n"
            "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
            "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))).. (-17.50)\n";
    std::stringstream istream{input};
    seqan3::structure_file_input fin{istream, seqan3::format_vienna{}, seqan3::fields<seqan3::field::structured_seq>{}};
    EXPECT_THROW(fin.begin(), seqan3::parse_error);
}

TEST_F(read_fail, structure_too_short_structured_seq)
{
    input = "> S.cerevisiae_tRNA-PHE M10740/1-73\n"
            "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
            "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))) (-17.50)\n";
    std::stringstream istream{input};
    seqan3::structure_file_input fin{istream, seqan3::format_vienna{}, seqan3::fields<seqan3::field::structured_seq>{}};
    EXPECT_THROW(fin.begin(), seqan3::parse_error);
}

TEST_F(read_fail, wrong_char)
{
    input = "> S.cerevisiae_tRNA-PHE M10740/1-73\n"
            "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUICUGUGUUCGAUCCACAGAAUUCGCA\n"
            "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)\n";
    std::stringstream istream{input};
    seqan3::structure_file_input fin{istream, seqan3::format_vienna{}};
    EXPECT_THROW(fin.begin(), seqan3::parse_error);
}

// ----------------------------------------------------------------------------
// writing
// ----------------------------------------------------------------------------

struct write : public ::testing::Test
{
    std::vector<std::string> id{{"S.cerevisiae_tRNA-PHE M10740/1-73"}, {"example 2"}};

    std::vector<seqan3::rna5_vector> const seq{
        {"GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA"_rna5},
        {"UUGGAGUACACAACCUGUACACUCUUUC"_rna5}};

    std::vector<std::vector<seqan3::wuss<51>>> const structure{
        {"(((((((..((((........)))).((((.........)))).....(((((.......))))))))))))."_wuss51},
        {"..(((((..(((...)))..)))))..."_wuss51}};

    std::vector<float> const energy{-17.5f, -3.71f};

    std::ostringstream ostream;
};

TEST_F(write, standard)
{
    seqan3::structure_file_output fout{
        ostream,
        seqan3::format_vienna{},
        seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::structure, seqan3::field::energy>{}};
    for (size_t i = 0ul; i < seq.size(); ++i)
        fout.emplace_back(seq[i], id[i], structure[i], energy[i]);
    //EXPECT_NO_THROW(format.write(ostream, options, seq[i], id[i], ig, structure[i], energy[i], ig, ig, ig, ig));

    std::string const expected_content{
        "> S.cerevisiae_tRNA-PHE M10740/1-73\n"
        "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
        "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.500000)\n"
        "> example 2\n"
        "UUGGAGUACACAACCUGUACACUCUUUC\n"
        "..(((((..(((...)))..)))))... (-3.710000)\n"};
    ostream.flush();
    EXPECT_EQ(ostream.str(), expected_content);
}

TEST_F(write, option_precision)
{
    seqan3::structure_file_output fout{
        ostream,
        seqan3::format_vienna{},
        seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::structure, seqan3::field::energy>{}};
    fout.options.precision = 2; // we want 2 digits for energy
    for (size_t i = 0ul; i < seq.size(); ++i)
        fout.emplace_back(seq[i], id[i], structure[i], energy[i]);

    std::string const expected_content{
        "> S.cerevisiae_tRNA-PHE M10740/1-73\n"
        "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
        "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)\n"
        "> example 2\n"
        "UUGGAGUACACAACCUGUACACUCUUUC\n"
        "..(((((..(((...)))..)))))... (-3.71)\n"};
    ostream.flush();
    EXPECT_EQ(ostream.str(), expected_content);
}

TEST_F(write, option_add_carriage_return)
{
    seqan3::structure_file_output fout{
        ostream,
        seqan3::format_vienna{},
        seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::structure, seqan3::field::energy>{}};
    fout.options.add_carriage_return = true;
    for (size_t i = 0ul; i < seq.size(); ++i)
        fout.emplace_back(seq[i], id[i], structure[i], energy[i]);

    std::string const expected_content{
        "> S.cerevisiae_tRNA-PHE M10740/1-73\r\n"
        "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\r\n"
        "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.500000)\r\n"
        "> example 2\r\n"
        "UUGGAGUACACAACCUGUACACUCUUUC\r\n"
        "..(((((..(((...)))..)))))... (-3.710000)\r\n"};
    ostream.flush();
    EXPECT_EQ(ostream.str(), expected_content);
}

struct write_fields : public write
{};

TEST_F(write_fields, id_missing)
{
    seqan3::structure_file_output fout{
        ostream,
        seqan3::format_vienna{},
        seqan3::fields<seqan3::field::seq, seqan3::field::structure, seqan3::field::energy>{}};
    fout.options.precision = 2; // we want 2 digits for energy
    for (size_t i = 0ul; i < seq.size(); ++i)
        fout.emplace_back(seq[i], structure[i], energy[i]);

    std::string const expected_content{
        "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
        "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)\n"
        "UUGGAGUACACAACCUGUACACUCUUUC\n"
        "..(((((..(((...)))..)))))... (-3.71)\n"};
    ostream.flush();
    EXPECT_EQ(ostream.str(), expected_content);
}

TEST_F(write_fields, energy_missing)
{
    seqan3::structure_file_output fout{
        ostream,
        seqan3::format_vienna{},
        seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::structure>{}};
    for (size_t i = 0ul; i < seq.size(); ++i)
        fout.emplace_back(seq[i], id[i], structure[i]);

    std::string const expected_content{"> S.cerevisiae_tRNA-PHE M10740/1-73\n"
                                       "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
                                       "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))).\n"
                                       "> example 2\n"
                                       "UUGGAGUACACAACCUGUACACUCUUUC\n"
                                       "..(((((..(((...)))..)))))...\n"};
    ostream.flush();
    EXPECT_EQ(ostream.str(), expected_content);
}

TEST_F(write_fields, structure_missing)
{
    seqan3::structure_file_output fout{ostream,
                                       seqan3::format_vienna{},
                                       seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::energy>{}};
    EXPECT_THROW(fout.emplace_back(seq[0], id[0], energy[0]), std::logic_error);
}

TEST_F(write_fields, structure_and_energy_missing)
{
    seqan3::structure_file_output fout{ostream,
                                       seqan3::format_vienna{},
                                       seqan3::fields<seqan3::field::seq, seqan3::field::id>{}};
    for (size_t i = 0ul; i < seq.size(); ++i)
        fout.emplace_back(seq[i], id[i]);

    std::string const expected_content{"> S.cerevisiae_tRNA-PHE M10740/1-73\n"
                                       "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
                                       "> example 2\n"
                                       "UUGGAGUACACAACCUGUACACUCUUUC\n"};
    ostream.flush();
    EXPECT_EQ(ostream.str(), expected_content);
}

TEST_F(write_fields, seq_missing)
{
    seqan3::structure_file_output fout{
        ostream,
        seqan3::format_vienna{},
        seqan3::fields<seqan3::field::id, seqan3::field::structure, seqan3::field::energy>{}};
    EXPECT_THROW(fout.emplace_back(id[0], structure[0], energy[0]), std::logic_error);
}

TEST_F(write_fields, seq_empty)
{
    seqan3::structure_file_output fout{
        ostream,
        seqan3::format_vienna{},
        seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::structure, seqan3::field::energy>{}};
    EXPECT_THROW(fout.emplace_back(""_rna5, id[0], structure[0], energy[0]), std::runtime_error);
}

TEST_F(write_fields, only_seq)
{
    seqan3::structure_file_output fout{ostream, seqan3::format_vienna{}, seqan3::fields<seqan3::field::seq>{}};
    for (size_t i = 0ul; i < seq.size(); ++i)
        fout.emplace_back(seq[i]);

    std::string const expected_content{"GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
                                       "UUGGAGUACACAACCUGUACACUCUUUC\n"};
    ostream.flush();
    EXPECT_EQ(ostream.str(), expected_content);
}

TEST_F(write_fields, structured_seq)
{
    std::vector<std::vector<seqan3::structured_rna<seqan3::rna5, seqan3::wuss<51>>>> structured_seq;
    structured_seq.resize(seq.size());

    for (unsigned i = 0; i < seq.size(); ++i)
    {
        structured_seq[i].resize(seq[i].size());
        for (unsigned j = 0; j < seq[i].size(); ++j)
        {
            structured_seq[i][j] = seq[i][j];
            structured_seq[i][j] = structure[i][j];
        }
    }

    seqan3::structure_file_output fout{ostream,
                                       seqan3::format_vienna{},
                                       seqan3::fields<seqan3::field::structured_seq>{}};

    for (size_t i = 0ul; i < seq.size(); ++i)
        fout.emplace_back(structured_seq[i]);

    std::string const expected_content{"GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
                                       "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))).\n"
                                       "UUGGAGUACACAACCUGUACACUCUUUC\n"
                                       "..(((((..(((...)))..)))))...\n"};
    ostream.flush();
    EXPECT_EQ(ostream.str(), expected_content);
}
