// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <tuple>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/io/detail/record.hpp>
#include <seqan3/test/expect_same_type.hpp>
#include <seqan3/utility/views/zip.hpp>

TEST(detail, select_types_with_ids)
{
    using types = seqan3::type_list<std::string, seqan3::dna4_vector, std::vector<seqan3::phred42>>;
    using types_as_ids = seqan3::fields<seqan3::field::id, seqan3::field::seq, seqan3::field::qual>;
    using selected_ids = seqan3::fields<seqan3::field::qual, seqan3::field::id>;

    using selected_types = typename seqan3::detail::select_types_with_ids<types, types_as_ids, selected_ids>::type;

    EXPECT_TRUE((std::is_same_v<selected_types, seqan3::type_list<std::vector<seqan3::phred42>, std::string>>));
}

TEST(get_or_ignore, record)
{
    using types = seqan3::type_list<std::string, seqan3::dna4_vector>;
    using types_as_ids = seqan3::fields<seqan3::field::id, seqan3::field::seq>;
    using record_type = seqan3::record<types, types_as_ids>;
    record_type record{};

    EXPECT_SAME_TYPE(std::string &, decltype(seqan3::detail::get_or_ignore<seqan3::field::id>(record)));
    EXPECT_SAME_TYPE(seqan3::dna4_vector &, decltype(seqan3::detail::get_or_ignore<seqan3::field::seq>(record)));
    EXPECT_SAME_TYPE(decltype(std::ignore) &, decltype(seqan3::detail::get_or_ignore<seqan3::field::qual>(record)));

    EXPECT_SAME_TYPE(std::string const &,
                     decltype(seqan3::detail::get_or_ignore<seqan3::field::id>(std::as_const(record))));
    EXPECT_SAME_TYPE(seqan3::dna4_vector const &,
                     decltype(seqan3::detail::get_or_ignore<seqan3::field::seq>(std::as_const(record))));
    EXPECT_SAME_TYPE(decltype(std::ignore) const &,
                     decltype(seqan3::detail::get_or_ignore<seqan3::field::qual>(std::as_const(record))));
}

TEST(get_or_ignore, tuple)
{
    std::tuple<std::string, seqan3::dna4_vector> tuple{};

    EXPECT_SAME_TYPE(std::string &, decltype(seqan3::detail::get_or_ignore<0>(tuple)));
    EXPECT_SAME_TYPE(seqan3::dna4_vector &, decltype(seqan3::detail::get_or_ignore<1>(tuple)));
    EXPECT_SAME_TYPE(decltype(std::ignore) &, decltype(seqan3::detail::get_or_ignore<2>(tuple)));

    EXPECT_SAME_TYPE(std::string const &, decltype(seqan3::detail::get_or_ignore<0>(std::as_const(tuple))));
    EXPECT_SAME_TYPE(seqan3::dna4_vector const &, decltype(seqan3::detail::get_or_ignore<1>(std::as_const(tuple))));
    EXPECT_SAME_TYPE(decltype(std::ignore) const &, decltype(seqan3::detail::get_or_ignore<2>(std::as_const(tuple))));
}

TEST(get_or_ignore, zip_tuple)
{
    std::vector<std::string> ids{{}};
    std::vector<seqan3::dna4_vector> sequences{{}};

    auto id_sequence_zip = seqan3::views::zip(ids, sequences);
    auto tuple = *std::ranges::begin(id_sequence_zip);

    EXPECT_SAME_TYPE(std::string &, decltype(seqan3::detail::get_or_ignore<0>(tuple)));
    EXPECT_SAME_TYPE(seqan3::dna4_vector &, decltype(seqan3::detail::get_or_ignore<1>(tuple)));
    EXPECT_SAME_TYPE(decltype(std::ignore) &, decltype(seqan3::detail::get_or_ignore<2>(tuple)));

    EXPECT_SAME_TYPE(std::string const &, decltype(seqan3::detail::get_or_ignore<0>(std::as_const(tuple))));
    EXPECT_SAME_TYPE(seqan3::dna4_vector const &, decltype(seqan3::detail::get_or_ignore<1>(std::as_const(tuple))));
    EXPECT_SAME_TYPE(decltype(std::ignore) const &, decltype(seqan3::detail::get_or_ignore<2>(std::as_const(tuple))));
}

template <typename id_t>
struct custom_tuple : public std::tuple<id_t, seqan3::dna4_vector>
{
    using base_t = std::tuple<id_t, seqan3::dna4_vector>;
};

namespace std
{
template <typename id_t>
struct tuple_size<custom_tuple<id_t>> : public tuple_size<typename custom_tuple<id_t>::base_t>
{};

template <size_t index, typename id_t>
struct tuple_element<index, custom_tuple<id_t>> : public tuple_element<index, typename custom_tuple<id_t>::base_t>
{};
} // namespace std

TEST(get_or_ignore, custom_tuple)
{
    custom_tuple<std::string> tuple{};

    EXPECT_SAME_TYPE(std::string &, decltype(seqan3::detail::get_or_ignore<0>(tuple)));
    EXPECT_SAME_TYPE(seqan3::dna4_vector &, decltype(seqan3::detail::get_or_ignore<1>(tuple)));
    EXPECT_SAME_TYPE(decltype(std::ignore) &, decltype(seqan3::detail::get_or_ignore<2>(tuple)));

    EXPECT_SAME_TYPE(std::string const &, decltype(seqan3::detail::get_or_ignore<0>(std::as_const(tuple))));
    EXPECT_SAME_TYPE(seqan3::dna4_vector const &, decltype(seqan3::detail::get_or_ignore<1>(std::as_const(tuple))));
    EXPECT_SAME_TYPE(decltype(std::ignore) const &, decltype(seqan3::detail::get_or_ignore<2>(std::as_const(tuple))));
}
