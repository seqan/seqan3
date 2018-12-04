// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2018, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <gtest/gtest.h>

#include <seqan3/core/metafunction/template_inspection.hpp>

namespace seqan3
{

//!\cond DEV
/*!\brief Transforms any type list into a google type list, i.e. ::%testing::Types.
 * \ingroup core
 * \tparam  type_list  The type list.
 *
 * \details
 *
 * You should use this if you would have to redefine every type in a second type list just to test the same types in
 * gtest's `TYPED_TEST` and in one of your test cases. This reduces redundancy and errors when keeping the types
 * up-to-date and in sync.
 *
 * \attention
 * Google type lists can't easily be transformed by seqan3::detail::transfer_template_args_onto to a seqan3::type_list,
 * because they are defined by 50 template arguments which are defaulted to `::%testing::internal::None`. That means that
 * `::%testing::Types<int>` would yield `seqan3::type_list<int, ::%testing::internal::None, ::%testing::internal::None,
 * ..., ::%testing::internal::None>` while you would only expect `seqan3::type_list<int>`. Therefore, first define the
 * seqan3::type_list and use seqan3::as_testing_list to transform it to a google type list.
 *
 * ## Example
 *
 * \include test/snippet/test/testing_list.cpp
 */
template <typename type_list>
using as_testing_list = detail::transfer_template_args_onto_t<type_list, ::testing::Types>;
//!\endcond

} // namespace seqan3
