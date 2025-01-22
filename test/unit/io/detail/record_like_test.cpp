// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/io/detail/record_like.hpp>

TEST(record_like, record)
{
    using types = seqan3::type_list<std::string, std::string>;
    using types_as_ids = seqan3::fields<seqan3::field::id, seqan3::field::seq>;
    using record_type = seqan3::record<types, types_as_ids>;

    EXPECT_FALSE((seqan3::detail::record_like<types>));
    EXPECT_FALSE((seqan3::detail::record_like<types_as_ids>));
    EXPECT_TRUE((seqan3::detail::record_like<record_type>));
}

struct my_record :
    seqan3::record<seqan3::type_list<std::string, std::string>, seqan3::fields<seqan3::field::id, seqan3::field::seq>>
{};

namespace std
{
template <>
struct tuple_size<my_record> : tuple_size<typename my_record::base_type>
{};

template <size_t id>
struct tuple_element<id, my_record> : tuple_element<id, typename my_record::base_type>
{};
} // namespace std

TEST(record_like, custom_record)
{
    EXPECT_TRUE((seqan3::detail::record_like<my_record>));
}
