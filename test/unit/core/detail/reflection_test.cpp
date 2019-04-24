// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/core/detail/reflection.hpp>

using namespace std::literals;
using namespace seqan3;

// Some test struct to get the name for.
namespace foo
{
template <typename ...type>
struct bar
{};
} // namespace foo

// Some types to test if reflection works.
using reflection_types = ::testing::Types<char, char16_t, char32_t, short,
                                          short int, unsigned int, double, const char *,
                                          foo::bar<char>, foo::bar<foo::bar<char, double>>,
                                          foo::bar<foo::bar<char, short *>>>;

template <typename param_type>
class reflection : public ::testing::Test
{
    // The corresponding list of names that should be generated. Must have the same order as `reflection_types`.
    inline static const std::vector names{"char", "char16_t", "char32_t", "short int", "short int", "unsigned int", "double",
                                          "const char*", "foo::bar<char>", "foo::bar<foo::bar<char, double> >",
                                          "foo::bar<foo::bar<char, short int*> >"};

    // Helper function to obtain the index of the current type `param_type` within the list of `reflection_types`.
    template <typename ...rem_types>
    auto type_to_index(::testing::Types<param_type, rem_types...> const & /*tuple*/)
    {
        return 0;
    }

    template <typename first_type, typename ...rem_types>
    size_t type_to_index(::testing::Types<first_type, rem_types...> const & /*tuple*/)
    {
        return 1 + type_to_index(::testing::Types<rem_types...>{});
    }

public:

    // Public interface to query the expected name for the current test instance.
    std::string expected_name()
    {
        return names[type_to_index(reflection_types{})];
    }
};

TYPED_TEST_CASE(reflection, reflection_types);

TYPED_TEST(reflection, size)
{
    EXPECT_EQ(detail::get_display_name_size_v<TypeParam>, this->expected_name().size());
}

TYPED_TEST(reflection, name)
{
    EXPECT_EQ(detail::get_display_name_v<TypeParam>.string(), this->expected_name());
}

namespace seqan3::detail
{
template <typename type>
struct pretty_function
{
    static auto to_string()
    {
        return std::string{__PRETTY_FUNCTION__};
    }
};
} // namespace seqan3::detail

// NOTE: This does not test a seqan3 library feature, but the underlying magic
// of how seqan3::detail::get_display_name_v works. This is needed, because
// __PRETTY_FUNCTION__ has vendor specific output which is not standardised.
TEST(pretty_function, to_string)
{
    std::string int_name = seqan3::detail::pretty_function<int>::to_string();
#if defined(__clang__)
    EXPECT_EQ(int_name, "static auto seqan3::detail::pretty_function<int>::to_string() [type = int]");
#elif defined(__GNUC__)
    EXPECT_EQ(int_name, "static auto seqan3::detail::pretty_function<type>::to_string() [with type = int]");
#endif

    std::string foo_bar_char_name = seqan3::detail::pretty_function<foo::bar<char>>::to_string();
#if defined(__clang__)
    EXPECT_EQ(foo_bar_char_name, "static auto seqan3::detail::pretty_function<foo::bar<char> >::to_string() [type = foo::bar<char>]");
#elif defined(__GNUC__)
    EXPECT_EQ(foo_bar_char_name, "static auto seqan3::detail::pretty_function<type>::to_string() [with type = foo::bar<char>]");
#endif
}
