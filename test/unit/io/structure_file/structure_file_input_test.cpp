// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <iterator>
#include <fstream>
#include <sstream>

#include <gtest/gtest.h>

#include <range/v3/view/zip.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/view/map.hpp>

#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/io/structure_file/input.hpp>
#include <seqan3/range/view/convert.hpp>
#include <seqan3/std/iterator>
#include <seqan3/std/view/filter.hpp>
#include <seqan3/test/tmp_filename.hpp>

using namespace seqan3;

TEST(general, concepts)
{
    using it_t = typename structure_file_in<>::iterator;
    using sen_t = typename structure_file_in<>::sentinel;

    EXPECT_TRUE((std::InputIterator<it_t>));
    EXPECT_TRUE((std::Sentinel<sen_t, it_t>));
}

struct structure_file_input_class : public ::testing::Test
{
    using comp0 = structure_file_input_default_traits_rna;
    using comp1 = fields<field::SEQ, field::ID, field::STRUCTURE>;
    using comp2 = type_list<structure_file_format_vienna>;
    using comp3 = std::ifstream;

    test::tmp_filename create_file(char const * contents)
    {
        test::tmp_filename filename{"structure_file_input_constructor.dbn"};
        {
            std::ofstream filecreator{filename.get_path(), std::ios::out | std::ios::binary};
            filecreator << contents; // must contain at least one record
        }
        return std::move(filename);
    }
};

TEST_F(structure_file_input_class, concepts)
{
    using t = structure_file_in<>;
    EXPECT_TRUE((std::ranges::InputRange<t>));

    using ct = structure_file_in<> const;
    // not const-iterable
    EXPECT_FALSE((std::ranges::InputRange<ct>));
}

TEST_F(structure_file_input_class, construct_by_filename)
{
    /* just the filename */
    {
        test::tmp_filename filename = create_file("> ID\nACGU\n....\n");
        EXPECT_NO_THROW(structure_file_in<>{filename.get_path()});
    }

    // correct format check is done by tests of that format

    /* wrong extension */
    {
        test::tmp_filename filename{"structure_file_input_constructor.xyz"};
        std::ofstream filecreator{filename.get_path(), std::ios::out | std::ios::binary};
        EXPECT_THROW(structure_file_in<>{filename.get_path()}, unhandled_extension_error);
    }

    /* non-existent file*/
    {
        EXPECT_THROW(structure_file_in<>{"/dev/nonexistant/foobarOOO"}, file_open_error);
    }

    /* filename + fields */
    {
        test::tmp_filename filename = create_file("> ID\nACGU\n....\n");
        EXPECT_NO_THROW((structure_file_in<structure_file_input_default_traits_rna,
                                           fields<field::SEQ>,
                                           type_list<structure_file_format_vienna>,
                                           std::ifstream>{filename.get_path(), fields<field::SEQ>{}}));
    }
}

TEST_F(structure_file_input_class, construct_from_stream)
{
    /* stream + format_tag */
    EXPECT_NO_THROW((structure_file_in<structure_file_input_default_traits_rna,
                                       fields<field::SEQ, field::ID, field::STRUCTURE>,
                                       type_list<structure_file_format_vienna>,
                                       std::istringstream>{std::istringstream{"> ID\nACGU\n....\n"},
                                                           structure_file_format_vienna{}}));


    /* stream + format_tag + fields */
    EXPECT_NO_THROW((structure_file_in<structure_file_input_default_traits_rna,
                                       fields<field::SEQ, field::ID, field::STRUCTURE>,
                                       type_list<structure_file_format_vienna>,
                                       std::istringstream>{std::istringstream{"> ID\nACGU\n....\n"},
                                                           structure_file_format_vienna{},
                                                           fields<field::SEQ, field::ID, field::STRUCTURE>{}}));
}

TEST_F(structure_file_input_class, default_template_args)
{
    /* default template args */
    using t = structure_file_in<>;
    EXPECT_TRUE((std::is_same_v<typename t::traits_type,        comp0>));
    EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
    EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      comp2>));
    EXPECT_TRUE((std::is_same_v<typename t::stream_type,        comp3>));
}

TEST_F(structure_file_input_class, guided_filename_constructor)
{
    /* guided filename constructor */
    test::tmp_filename filename = create_file("> ID\nACGU\n....\n");
    structure_file_in fin{filename.get_path()};

    using t = decltype(fin);
    EXPECT_TRUE((std::is_same_v<typename t::traits_type, comp0>));
    EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
    EXPECT_TRUE((std::is_same_v<typename t::valid_formats, comp2>));
    EXPECT_TRUE((std::is_same_v<typename t::stream_type, comp3>));
}

TEST_F(structure_file_input_class, guided_filename_constructor_and_custom_fields)
{
    /* guided filename constructor + custom fields */
    test::tmp_filename filename = create_file("> ID\nACGU\n....\n");
    structure_file_in fin{filename.get_path(), fields<field::SEQ>{}};

    using t = decltype(fin);
    EXPECT_TRUE((std::is_same_v<typename t::traits_type,        comp0>));
    EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, fields<field::SEQ>>));
    EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      comp2>));
    EXPECT_TRUE((std::is_same_v<typename t::stream_type,        comp3>));
}

TEST_F(structure_file_input_class, guided_stream_constructor)
{
    /* guided stream constructor */
    structure_file_in fin{std::istringstream{"> ID\nACGU\n....\n"}, structure_file_format_vienna{}};

    using t = decltype(fin);
    EXPECT_TRUE((std::is_same_v<typename t::traits_type,        comp0>));
    EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
    EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      type_list<structure_file_format_vienna>>));
    EXPECT_TRUE((std::is_same_v<typename t::stream_type,        std::istringstream>));
}

TEST_F(structure_file_input_class, guided_stream_constructor_and_custom_fields)
{
    /* guided stream constructor + custom fields */
    structure_file_in fin{std::istringstream{"> ID\nACGU\n....\n"},
                          structure_file_format_vienna{},
                          fields<field::SEQ>{}};

    using t = decltype(fin);
    EXPECT_TRUE((std::is_same_v<typename t::traits_type,        comp0>));
    EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, fields<field::SEQ>>));
    EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      type_list<structure_file_format_vienna>>));
    EXPECT_TRUE((std::is_same_v<typename t::stream_type,        std::istringstream>));
}

TEST_F(structure_file_input_class, amino_acids_traits)
{
    test::tmp_filename filename = create_file("> ID\nACEW\nHHHH\n");
    structure_file_in<structure_file_input_default_traits_aa> fin{filename.get_path()};

    using t = decltype(fin);
    EXPECT_TRUE((std::is_same_v<typename t::traits_type,        structure_file_input_default_traits_aa>));
    EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
    EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      comp2>));
    EXPECT_TRUE((std::is_same_v<typename t::stream_type,        comp3>));
}

TEST_F(structure_file_input_class, modified_traits)
{
    test::tmp_filename filename = create_file("> ID\nACGU\n....\n");
    //! [structure_file_input_class mod_traits]
    struct my_traits : structure_file_input_default_traits_rna
    {
        using seq_alphabet = rna4; // instead of rna5
    };

    structure_file_in<my_traits> fin{filename.get_path()};
    //! [structure_file_input_class mod_traits]

    using t = decltype(fin);
    EXPECT_TRUE((std::is_same_v<typename t::traits_type,        my_traits>));
    EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
    EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      comp2>));
    EXPECT_TRUE((std::is_same_v<typename t::stream_type,        comp3>));
}

struct structure_file_input_read : public ::testing::Test
{
    size_t const num_records = 2ul;

    std::string const input
    {
        ">S.cerevisiae_tRNA-PHE M10740/1-73\n"
        "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
        "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)\n"
        "> example\n"
        "UUGGAGUACACAACCUGUACACUCUUUC\n"
        "..(((((..(((...)))..)))))... (-3.71)\n"
    };

    rna5_vector const seq_comp[2]
    {
        "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA"_rna5,
        "UUGGAGUACACAACCUGUACACUCUUUC"_rna5
    };

    std::string const id_comp[2]
    {
        "S.cerevisiae_tRNA-PHE M10740/1-73",
        "example"
    };

    double const energy_comp[2]
    {
        -17.5, -3.71
    };

    std::vector<wuss51> const structure_comp[2]
    {
        "(((((((..((((........)))).((((.........)))).....(((((.......))))))))))))."_wuss51,
        "..(((((..(((...)))..)))))..."_wuss51
    };

    std::vector<uint8_t> const interaction_comp[2]
    {
        {
            71, 70, 69, 68, 67, 66, 65, 24, 23, 22, 21, 12, 11, 10,  9, 42, 41, 40, 39, 29,
            28, 27, 26, 64, 63, 62, 61, 60, 52, 51, 50, 49, 48,  6,  5,  4,  3,  2,  1,  0
        },
        {
            24, 23, 22, 21, 20, 17, 16, 15, 11, 10,  9,  6,  5,  4,  3,  2
        }
    };

    template <typename bpp_type>
    void bpp_test(bpp_type & bpp, std::vector<uint8_t> const & bpp_comp)
    {
        size_t idx = 0ul;
        auto interactions = bpp | ranges::view::remove_if([] (auto & set) { return set.size() != 1; });
        for (auto & elem : interactions)
        {
            EXPECT_EQ(elem.begin()->second, bpp_comp[idx++]);
        }
        EXPECT_EQ(idx, bpp_comp.size());
    }
};

TEST_F(structure_file_input_read, record_general)
{
    /* record based reading */
    structure_file_in fin{std::istringstream{input}, structure_file_format_vienna{}};

    size_t counter = 0ul;
    for (auto & rec : fin)
    {
        EXPECT_TRUE((std::ranges::equal(get<field::SEQ>(rec), seq_comp[counter])));
        EXPECT_TRUE((std::ranges::equal(get<field::ID>(rec), id_comp[counter])));
        EXPECT_TRUE((std::ranges::equal(get<field::STRUCTURE>(rec), structure_comp[counter])));
        ++counter;
    }
    EXPECT_EQ(counter, num_records);
}

TEST_F(structure_file_input_read, record_struct_bind)
{
    /* record based reading */
    structure_file_in fin{std::istringstream{input},
                          structure_file_format_vienna{},
                          fields<field::SEQ, field::ID, field::BPP, field::STRUCTURE, field::ENERGY>{}};

    size_t counter = 0ul;
    for (auto & [ sequence, id, bpp, structure, energy ] : fin)
    {
        EXPECT_TRUE((std::ranges::equal(sequence, seq_comp[counter])));
        EXPECT_TRUE((std::ranges::equal(id, id_comp[counter])));
        EXPECT_TRUE((std::ranges::equal(structure, structure_comp[counter])));
        EXPECT_DOUBLE_EQ(energy.value(), energy_comp[counter]);
        bpp_test(bpp, interaction_comp[counter]);
        ++counter;
    }
    EXPECT_EQ(counter, num_records);
}

TEST_F(structure_file_input_read, record_custom_fields)
{
    /* record based reading */
    structure_file_in fin{std::istringstream{input},
                         structure_file_format_vienna{},
                         fields<field::ID, field::STRUCTURED_SEQ>{}};

    size_t counter = 0ul;
    for (auto & [ id, seq_structure ] : fin)
    {
        EXPECT_TRUE((std::ranges::equal(id, id_comp[counter])));
        EXPECT_TRUE((std::ranges::equal(seq_structure | view::convert<rna5>, seq_comp[counter])));
        EXPECT_TRUE((std::ranges::equal(seq_structure | view::convert<wuss51>, structure_comp[counter])));
        ++counter;
    }
    EXPECT_EQ(counter, num_records);
}

TEST_F(structure_file_input_read, record_file_view)
{
    structure_file_in fin{std::istringstream{input}, structure_file_format_vienna{},
                          fields<field::SEQ, field::ID, field::BPP, field::STRUCTURE, field::ENERGY>{}};

    auto minimum_length_filter = ranges::view::filter([] (auto const & rec)
    {
        return size(get<field::SEQ>(rec)) >= 5;
    });

    size_t counter = 0ul; // the first record will be filtered out
    for (auto & rec : fin | minimum_length_filter)
    {
        EXPECT_TRUE((std::ranges::equal(get<field::SEQ>(rec), seq_comp[counter])));
        EXPECT_TRUE((std::ranges::equal(get<field::ID>(rec),  id_comp[counter])));
        bpp_test(get<field::BPP>(rec), interaction_comp[counter]);
        EXPECT_TRUE((std::ranges::equal(get<field::STRUCTURE>(rec), structure_comp[counter])));
        EXPECT_DOUBLE_EQ(get<field::ENERGY>(rec).value(), energy_comp[counter]);
        ++counter;
    }
    EXPECT_EQ(counter, num_records);
}

TEST_F(structure_file_input_read, column_general)
{
    structure_file_in fin{std::istringstream{input}, structure_file_format_vienna{},
                          fields<field::SEQ, field::ID, field::BPP, field::STRUCTURE, field::ENERGY>{}};

    auto & seqs  = get<field::SEQ>(fin);                                    // by field
    auto & ids   = get<1>(fin);                                             // by index
    auto & bpps = get<field::BPP>(fin);
    auto & struc = get<typename decltype(fin)::structure_column_type>(fin); // by type
    auto & energies = get<field::ENERGY>(fin);

    ASSERT_EQ(seqs.size(), num_records);
    ASSERT_EQ(ids.size(), num_records);
    ASSERT_EQ(bpps.size(), num_records);
    ASSERT_EQ(struc.size(), num_records);
    ASSERT_EQ(energies.size(), num_records);

    for (size_t idx = 0ul; idx < num_records; ++idx)
    {
        EXPECT_TRUE((std::ranges::equal(seqs[idx], seq_comp[idx])));
        EXPECT_TRUE((std::ranges::equal(ids[idx], id_comp[idx])));
        bpp_test(bpps[idx], interaction_comp[idx]);
        EXPECT_TRUE((std::ranges::equal(struc[idx], structure_comp[idx])));
        EXPECT_DOUBLE_EQ(energies[idx].value(), energy_comp[idx]);
    }
}

TEST_F(structure_file_input_read, column_temporary)
{
    structure_file_in{std::istringstream{input}, structure_file_format_vienna{}};

    auto seqs = get<field::SEQ>(structure_file_in{std::istringstream{input}, structure_file_format_vienna{}});

    ASSERT_EQ(seqs.size(), num_records);

    for (size_t idx = 0ul; idx < num_records; ++idx)
    {
        EXPECT_TRUE((std::ranges::equal(seqs[idx], seq_comp[idx])));
    }
}

TEST_F(structure_file_input_read, column_decomposed)
{
    structure_file_in fin{std::istringstream{input}, structure_file_format_vienna{},
                          fields<field::SEQ, field::ID, field::STRUCTURE, field::ENERGY, field::BPP>{}};

    auto & [ seqs, ids, struc, energies, bpps ] = fin;

    ASSERT_EQ(seqs.size(), num_records);
    ASSERT_EQ(ids.size(), num_records);
    ASSERT_EQ(struc.size(), num_records);
    ASSERT_EQ(energies.size(), num_records);
    ASSERT_EQ(bpps.size(), num_records);

    for (size_t idx = 0ul; idx < num_records; ++idx)
    {
        EXPECT_TRUE((std::ranges::equal(seqs[idx], seq_comp[idx])));
        EXPECT_TRUE((std::ranges::equal(ids[idx], id_comp[idx])));
        EXPECT_TRUE((std::ranges::equal(struc[idx], structure_comp[idx])));
        EXPECT_DOUBLE_EQ(energies[idx].value(), energy_comp[idx]);
        bpp_test(bpps[idx], interaction_comp[idx]);
    }
}

TEST_F(structure_file_input_read, column_decomposed_temporary)
{
    auto && [ seqs, ids, struc, energies, bpps ] = structure_file_in{std::istringstream{input},
                                                                     structure_file_format_vienna{},
                                                                     fields<field::SEQ,
                                                                            field::ID,
                                                                            field::STRUCTURE,
                                                                            field::ENERGY,
                                                                            field::BPP>{}};

    ASSERT_EQ(seqs.size(), num_records);
    ASSERT_EQ(ids.size(), num_records);
    ASSERT_EQ(struc.size(), num_records);
    ASSERT_EQ(energies.size(), num_records);
    ASSERT_EQ(bpps.size(), num_records);

    for (size_t idx = 0ul; idx < num_records; ++idx)
    {
        EXPECT_TRUE((std::ranges::equal(seqs[idx], seq_comp[idx])));
        EXPECT_TRUE((std::ranges::equal(ids[idx], id_comp[idx])));
        EXPECT_TRUE((std::ranges::equal(struc[idx], structure_comp[idx])));
        EXPECT_DOUBLE_EQ(energies[idx].value(), energy_comp[idx]);
        bpp_test(bpps[idx], interaction_comp[idx]);
    }
}
