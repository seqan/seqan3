// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================

/*!\file
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 * \brief Tests for the Vienna file format.
 */

#include <optional>
#include <sstream>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include <range/v3/algorithm/equal.hpp>

#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/alphabet/structure/dot_bracket3.hpp>
#include <seqan3/alphabet/structure/structured_rna.hpp>
#include <seqan3/io/structure_file/format_vienna.hpp>
#include <seqan3/io/structure_file/input_format_concept.hpp>
#include <seqan3/io/structure_file/output_format_concept.hpp>
#include <seqan3/range/view/convert.hpp>

using namespace seqan3;

// ----------------------------------------------------------------------------
// general
// ----------------------------------------------------------------------------

TEST(general, concepts)
{
    EXPECT_TRUE((structure_file_input_format_concept<structure_file_format_vienna>));
    EXPECT_TRUE((structure_file_output_format_concept<structure_file_format_vienna>));
}

// ----------------------------------------------------------------------------
// reading
// ----------------------------------------------------------------------------

struct read : public ::testing::Test
{
    std::string input
    {
        "> S.cerevisiae_tRNA-PHE M10740/1-73\n"
        "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
        "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)\n"
        "> example 2\n"
        "UUGGAGUACACAACCUGUACACUCUUUC\n"
        "..(((((..(((...)))..)))))... (-3.71)\n"
    };

    std::vector<std::string> expected_id
    {
        { "S.cerevisiae_tRNA-PHE M10740/1-73" },
        { "example 2" }
    };

    std::vector<rna4_vector> const expected_seq
    {
        { "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA"_rna4 },
        { "UUGGAGUACACAACCUGUACACUCUUUC"_rna4 }
    };

    std::vector<std::vector<dot_bracket3>> const expected_structure
    {
        { "(((((((..((((........)))).((((.........)))).....(((((.......))))))))))))."_db3 },
        { "..(((((..(((...)))..)))))..."_db3 }
    };

    std::vector<double> const expected_energy
    {
        -17.5, -3.71
    };

    std::vector<std::vector<uint8_t>> const expected_interactions
    {
        {
            71, 70, 69, 68, 67, 66, 65, 24, 23, 22, 21, 12, 11, 10,  9, 42, 41, 40, 39, 29,
            28, 27, 26, 64, 63, 62, 61, 60, 52, 51, 50, 49, 48,  6,  5,  4,  3,  2,  1,  0
        },
        {
            24, 23, 22, 21, 20, 17, 16, 15, 11, 10,  9,  6,  5,  4,  3,  2
        }
    };

    structure_file_format_vienna format;

    structure_file_input_options<rna4, false> options;

    static constexpr auto ig = std::ignore; // shortcut

    bool check_seq = true;
    bool check_id = true;
    bool check_structure = true;
    bool check_energy= true;

    rna4_vector seq;
    std::string id;
    std::vector<std::set<std::pair<float, uint8_t>>> bpp;
    std::vector<dot_bracket3> structure;
    std::vector<structured_rna<rna4, dot_bracket3>> structured_seq;
    std::optional<double> energy;

    void bpp_test(std::vector<uint8_t> const & bpp_comp)
    {
        size_t cnt = 0ul;
        auto interaction_sets = bpp | ranges::view::remove_if([] (auto & set) { return set.size() != 1; });
        for (auto & iset : interaction_sets)
        {
            EXPECT_EQ(iset.size(), 1);
            EXPECT_EQ(iset.begin()->second, bpp_comp[cnt++]);
            EXPECT_FLOAT_EQ(iset.begin()->first, 1.0f);
        }
        EXPECT_EQ(cnt, bpp_comp.size());
    }

    void do_read_test(std::string const & input)
    {
        std::stringstream istream{input};

        for (size_t idx = 0ul; idx < expected_seq.size(); ++idx)
        {
            id.clear();
            seq.clear();
            structure.clear();
            bpp.clear();
            energy = std::nullopt;

            EXPECT_NO_THROW(( format.read(istream, options, seq, id, bpp, structure, energy, ig, ig, ig, ig) ));

            EXPECT_EQ(check_energy, energy.has_value());
            if (check_seq)       { EXPECT_TRUE(( std::ranges::equal(seq, expected_seq[idx]) )); }
            if (check_id)        { EXPECT_TRUE(( std::ranges::equal(id, expected_id[idx]) )); }
            if (check_structure) { EXPECT_TRUE(( std::ranges::equal(structure, expected_structure[idx]) )); }
            if (check_structure) { bpp_test(expected_interactions[idx]); }
            if (check_energy)    { EXPECT_DOUBLE_EQ(*energy, expected_energy[idx]); }
        }
    }
};

TEST_F(read, standard)
{
    do_read_test(input);
}

TEST_F(read, newline_before_eof)
{
    input =
    {
        "> S.cerevisiae_tRNA-PHE M10740/1-73\n"
        "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
        "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)\n"
        "> example 2\n"
        "UUGGAGUACACAACCUGUACACUCUUUC\n"
        "..(((((..(((...)))..)))))... (-3.71)"
    };
    do_read_test(input);
}

TEST_F(read, whitespace_in_seq)
{
    input =
    {
        "> S.cerevisiae_tRNA-PHE M10740/1-73\n"
        "GCGGAUUUAGCU CAGUUGGGAGAGCGCCAGACU GAAGAUUUGGAGGUC CUGUGUUCGAUCCACA   GAAUU CGCA\n"
        "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)\n"
        "> example 2\n"
        "UUGGAGUAC   ACAACCUGUACAC UCUUUC \n"
        "..(((((..(((...)))..)))))... (-3.71)\n"
    };
    do_read_test(input);
}

TEST_F(read, no_energies)
{
    check_energy = false;
    input =
    {
        "> S.cerevisiae_tRNA-PHE M10740/1-73\n"
        "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
        "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))).\n"
        "> example 2\n"
        "UUGGAGUACACAACCUGUACACUCUUUC\n"
        "..(((((..(((...)))..)))))...\n"
    };
    do_read_test(input);
}

TEST_F(read, no_ids)
{
    check_id = false;
    input =
    {
        "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
        "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)\n"
        "UUGGAGUACACAACCUGUACACUCUUUC\n"
        "..(((((..(((...)))..)))))... (-3.71)\n"
    };
    do_read_test(input);
}

TEST_F(read, spaces_and_carriage_return)
{
    check_id = false;
    input =
    {
        "> S.cerevisiae_tRNA-PHE M10740/1-73\r\n"
        "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\r\n"
        "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.5)\r\n"
        ">example 2\r\n"
        "UUGGAGUA CACAACCUGUACA  CUCU UUC \r\n"
        "..(((((..(((...)))..)))))...     ( -3.71 )\r\n"
    };
    do_read_test(input);
}

TEST_F(read, options_truncate_ids)
{
    input =
    {
        "> S.cerevisiae_tRNA-PHE M10740/1-73\n"
        "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
        "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)\n"
        "> example 2\n"
        "UUGGAGUACACAACCUGUACACUCUUUC\n"
        "..(((((..(((...)))..)))))... (-3.71)\n"
    };

    options.truncate_ids = true;
    expected_id = {"S.cerevisiae_tRNA-PHE", "example"};
    do_read_test(input);
}

struct read_fields : public read {};

TEST_F(read_fields, only_seq)
{
    std::stringstream istream{input};
    for (size_t idx = 0ul; idx < expected_seq.size(); ++idx)
    {
        EXPECT_NO_THROW(( format.read(istream, options, seq, ig, ig, ig, ig, ig, ig, ig, ig) ));

        EXPECT_TRUE(std::ranges::equal(seq, expected_seq[idx]));
        seq.clear();
    }
}

TEST_F(read_fields, only_id)
{
    std::stringstream istream{input};
    for (size_t idx = 0ul; idx < expected_seq.size(); ++idx)
    {
        EXPECT_NO_THROW(( format.read(istream, options, ig, id, ig, ig, ig, ig, ig, ig, ig) ));

        EXPECT_TRUE(std::ranges::equal(id, expected_id[idx]));
        id.clear();
    }
}

TEST_F(read_fields, only_structure)
{
    std::stringstream istream{input};
    for (size_t idx = 0ul; idx < expected_seq.size(); ++idx)
    {
        EXPECT_NO_THROW(( format.read(istream, options, ig, ig, ig, structure, ig, ig, ig, ig, ig) ));

        EXPECT_TRUE(std::ranges::equal(structure, expected_structure[idx]));
        structure.clear();
    }
}

TEST_F(read_fields, only_energy)
{
    std::stringstream istream{input};
    for (size_t idx = 0ul; idx < expected_seq.size(); ++idx)
    {
        EXPECT_NO_THROW(( format.read(istream, options, ig, ig, ig, ig, energy, ig, ig, ig, ig) ));

        EXPECT_TRUE(energy);
        EXPECT_DOUBLE_EQ(*energy, expected_energy[idx]);
        energy = std::nullopt;
    }
}

TEST_F(read_fields, structured_seq)
{
    structure_file_input_options<rna4, true> opt{};
    std::stringstream istream{input};
    for (size_t idx = 0ul; idx < expected_seq.size(); ++idx)
    {
        EXPECT_NO_THROW(( format.read(istream, opt, structured_seq, ig, ig, structured_seq, ig, ig, ig, ig, ig) ));

        EXPECT_TRUE(std::ranges::equal(structured_seq | view::convert<rna4>, expected_seq[idx]));
        EXPECT_TRUE(std::ranges::equal(structured_seq | view::convert<dot_bracket3>, expected_structure[idx]));
        structured_seq.clear();
    }
}

TEST_F(read_fields, only_bpp)
{
    std::stringstream istream{input};
    for (size_t idx = 0ul; idx < expected_seq.size(); ++idx)
    {
        EXPECT_NO_THROW(( format.read(istream, options, ig, ig, bpp, ig, ig, ig, ig, ig, ig) ));

        bpp_test(expected_interactions[idx]);
        bpp.clear();
    }
}

struct read_fail : public read {};

TEST_F(read_fail, wrong_id)
{
    input[0] = '#'; // invalid character for ID line
    std::stringstream istream{input};
    EXPECT_THROW( format.read(istream, options, ig, id, ig, ig, ig, ig, ig, ig, ig), parse_error );
}

TEST_F(read_fail, missing_seq)
{
    input =
    {
        "> S.cerevisiae_tRNA-PHE M10740/1-73\n"
        "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)\n"
    };
    std::stringstream istream{input};
    EXPECT_THROW( format.read(istream, options, seq, id, bpp, structure, energy, ig, ig, ig, ig), parse_error );
}

TEST_F(read_fail, missing_structure)
{
    input =
    {
        "> S.cerevisiae_tRNA-PHE M10740/1-73\n"
        "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
        "> example 2\n"
        "UUGGAGUACACAACCUGUACACUCUUUC\n"
    };
    std::stringstream istream{input};
    EXPECT_THROW( format.read(istream, options, seq, id, bpp, structure, energy, ig, ig, ig, ig), parse_error );
}

TEST_F(read_fail, missing_structure_and_id)
{
    input =
    {
        "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
        "UUGGAGUACACAACCUGUACACUCUUUC\n"
    };
    std::stringstream istream{input};
    EXPECT_THROW( format.read(istream, options, seq, id, bpp, structure, energy, ig, ig, ig, ig), parse_error );
}

TEST_F(read_fail, empty)
{
    input =
    {
        ""
    };
    std::stringstream istream{input};
    EXPECT_THROW( format.read(istream, options, seq, id, bpp, structure, energy, ig, ig, ig, ig), parse_error );
}

TEST_F(read_fail, structure_too_long)
{
    input =
    {
        "> S.cerevisiae_tRNA-PHE M10740/1-73\n"
        "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
        "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))).. (-17.50)\n"
    };
    std::stringstream istream{input};
    EXPECT_THROW( format.read(istream, options, seq, id, bpp, structure, energy, ig, ig, ig, ig),
                  parse_error );
}

TEST_F(read_fail, structure_too_short)
{
    input =
    {
        "> S.cerevisiae_tRNA-PHE M10740/1-73\n"
        "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
        "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))) (-17.50)\n"
    };
    std::stringstream istream{input};
    EXPECT_THROW( format.read(istream, options, seq, id, bpp, structure, energy, ig, ig, ig, ig),
                  parse_error );
}

TEST_F(read_fail, structure_too_long_structured_seq)
{
    input =
    {
        "> S.cerevisiae_tRNA-PHE M10740/1-73\n"
        "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
        "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))).. (-17.50)\n"
    };
    std::stringstream istream{input};
    EXPECT_THROW( format.read(istream, options, structured_seq, id, bpp, structured_seq, energy, ig, ig, ig, ig),
                  parse_error );
}

TEST_F(read_fail, structure_too_short_structured_seq)
{
    input =
    {
        "> S.cerevisiae_tRNA-PHE M10740/1-73\n"
        "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
        "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))) (-17.50)\n"
    };
    std::stringstream istream{input};
    EXPECT_THROW( format.read(istream, options, structured_seq, id, bpp, structured_seq, energy, ig, ig, ig, ig),
                  parse_error );
}

TEST_F(read_fail, wrong_char)
{
    input =
    {
        "> S.cerevisiae_tRNA-PHE M10740/1-73\n"
        "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUICUGUGUUCGAUCCACAGAAUUCGCA\n"
        "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)\n"
    };
    std::stringstream istream{input};
    EXPECT_THROW( format.read(istream, options, seq, id, bpp, structure, energy, ig, ig, ig, ig), parse_error );
}

// ----------------------------------------------------------------------------
// writing
// ----------------------------------------------------------------------------

struct write : public ::testing::Test
{
    std::vector<std::string> id
    {
        { "S.cerevisiae_tRNA-PHE M10740/1-73" },
        { "example 2" }
    };

    std::vector<rna4_vector> const seq
    {
        { "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA"_rna4 },
        { "UUGGAGUACACAACCUGUACACUCUUUC"_rna4 }
    };

    std::vector<std::vector<dot_bracket3>> const structure
    {
        { "(((((((..((((........)))).((((.........)))).....(((((.......))))))))))))."_db3 },
        { "..(((((..(((...)))..)))))..."_db3 }
    };

    std::vector<float> const energy
    {
        -17.5f, -3.71f
    };

    structure_file_format_vienna format;

    structure_file_output_options options;

    static constexpr auto ig = std::ignore; // shortcut

    std::ostringstream ostream;
};

TEST_F(write, standard)
{
    for (size_t i = 0ul; i < seq.size(); ++i)
        EXPECT_NO_THROW(format.write(ostream, options, seq[i], id[i], ig, structure[i], energy[i], ig, ig, ig, ig));

    std::string const expected_content
    {
        "> S.cerevisiae_tRNA-PHE M10740/1-73\n"
        "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
        "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.500000)\n"
        "> example 2\n"
        "UUGGAGUACACAACCUGUACACUCUUUC\n"
        "..(((((..(((...)))..)))))... (-3.710000)\n"
    };
    ostream.flush();
    EXPECT_EQ(ostream.str(), expected_content);
}

TEST_F(write, option_precision)
{
    options.precision = 2; // we want 2 digits for energy
    for (size_t i = 0ul; i < seq.size(); ++i)
        EXPECT_NO_THROW(format.write(ostream, options, seq[i], id[i], ig, structure[i], energy[i], ig, ig, ig, ig));

    std::string const expected_content
    {
        "> S.cerevisiae_tRNA-PHE M10740/1-73\n"
        "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
        "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)\n"
        "> example 2\n"
        "UUGGAGUACACAACCUGUACACUCUUUC\n"
        "..(((((..(((...)))..)))))... (-3.71)\n"
    };
    ostream.flush();
    EXPECT_EQ(ostream.str(), expected_content);
}

TEST_F(write, option_add_carriage_return)
{
    options.add_carriage_return = true;
    for (size_t i = 0ul; i < seq.size(); ++i)
        EXPECT_NO_THROW(format.write(ostream, options, seq[i], id[i], ig, structure[i], energy[i], ig, ig, ig, ig));

    std::string const expected_content
    {
        "> S.cerevisiae_tRNA-PHE M10740/1-73\r\n"
        "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\r\n"
        "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.500000)\r\n"
        "> example 2\r\n"
        "UUGGAGUACACAACCUGUACACUCUUUC\r\n"
        "..(((((..(((...)))..)))))... (-3.710000)\r\n"
    };
    ostream.flush();
    EXPECT_EQ(ostream.str(), expected_content);
}

struct write_fields : public write {};

TEST_F(write_fields, id_missing)
{
    options.precision = 2; // we want 2 digits for energy
    for (size_t i = 0ul; i < seq.size(); ++i)
        EXPECT_NO_THROW(format.write(ostream, options, seq[i], ig, ig, structure[i], energy[i], ig, ig, ig, ig));

    std::string const expected_content
    {
        "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
        "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)\n"
        "UUGGAGUACACAACCUGUACACUCUUUC\n"
        "..(((((..(((...)))..)))))... (-3.71)\n"
    };
    ostream.flush();
    EXPECT_EQ(ostream.str(), expected_content);
}

TEST_F(write_fields, energy_missing)
{
    for (size_t i = 0ul; i < seq.size(); ++i)
        EXPECT_NO_THROW(format.write(ostream, options, seq[i], id[i], ig, structure[i], ig, ig, ig, ig, ig));

    std::string const expected_content
    {
        "> S.cerevisiae_tRNA-PHE M10740/1-73\n"
        "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
        "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))).\n"
        "> example 2\n"
        "UUGGAGUACACAACCUGUACACUCUUUC\n"
        "..(((((..(((...)))..)))))...\n"
    };
    ostream.flush();
    EXPECT_EQ(ostream.str(), expected_content);
}

TEST_F(write_fields, structure_missing)
{
    EXPECT_THROW(format.write(ostream, options, seq[0], id[0], ig, ig, energy[0], ig, ig, ig, ig), std::logic_error);
}

TEST_F(write_fields, structure_and_energy_missing)
{
    for (size_t i = 0ul; i < seq.size(); ++i)
        EXPECT_NO_THROW(format.write(ostream, options, seq[i], id[i], ig, ig, ig, ig, ig, ig, ig));

    std::string const expected_content
    {
        "> S.cerevisiae_tRNA-PHE M10740/1-73\n"
        "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
        "> example 2\n"
        "UUGGAGUACACAACCUGUACACUCUUUC\n"
    };
    ostream.flush();
    EXPECT_EQ(ostream.str(), expected_content);
}

TEST_F(write_fields, seq_missing)
{
    EXPECT_THROW(format.write(ostream, options, ig, id[0], ig, structure[0], energy[0], ig, ig, ig, ig),
                 std::logic_error);
}

TEST_F(write_fields, seq_empty)
{
    EXPECT_THROW(format.write(ostream, options, ""_rna4, id[0], ig, structure[0], energy[0], ig, ig, ig, ig),
                 std::runtime_error);
}

TEST_F(write_fields, only_seq)
{
    for (size_t i = 0ul; i < seq.size(); ++i)
        EXPECT_NO_THROW(format.write(ostream, options, seq[i], ig, ig, ig, ig, ig, ig, ig, ig));

    std::string const expected_content
    {
        "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
        "UUGGAGUACACAACCUGUACACUCUUUC\n"
    };
    ostream.flush();
    EXPECT_EQ(ostream.str(), expected_content);
}

TEST_F(write_fields, structured_seq)
{
    std::vector<std::vector<structured_rna<rna4, dot_bracket3>>> structured_seq;
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

    for (size_t i = 0ul; i < seq.size(); ++i)
    {
        EXPECT_NO_THROW(format.write(ostream,
                                     options,
                                     structured_seq[i] | view::convert<rna4>,
                                     ig,
                                     ig,
                                     structured_seq[i] | view::convert<dot_bracket3>,
                                     ig,
                                     ig,
                                     ig,
                                     ig,
                                     ig));
    }

    std::string const expected_content
    {
        "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
        "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))).\n"
        "UUGGAGUACACAACCUGUACACUCUUUC\n"
        "..(((((..(((...)))..)))))...\n"
    };
    ostream.flush();
    EXPECT_EQ(ostream.str(), expected_content);
}
