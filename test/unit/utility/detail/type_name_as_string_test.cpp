// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/utility/detail/type_name_as_string.hpp>
#include <seqan3/utility/type_list/traits.hpp>
#include <seqan3/utility/type_list/type_list.hpp>

// Some test namespace to check if namespace information are preserved within the naming.
namespace foo
{
template <typename... type>
struct bar
{};
} // namespace foo

// Some types to test if type inspection works as expected.
// Note that the returned name might differ between compiler vendors and thus must be adapted accordingly
// in case this tests fails for those vendors.
using reflection_types = ::testing::Types<char,
                                          char16_t const,
                                          char32_t &,
                                          short *,
                                          double const * const,
                                          foo::bar<char> const &,
                                          foo::bar<foo::bar<char, double>>>;

// Helper type list to use some traits functions on type lists.
using as_type_list_t = seqan3::detail::transfer_template_args_onto_t<reflection_types, seqan3::type_list>;

template <typename param_type>
class type_inspection : public ::testing::Test
{
    // The corresponding list of names that should be generated. Must have the same order as `reflection_types`.
    static inline std::vector const names{"char",
                                          "char16_t const",
                                          "char32_t &",
                                          "short*",
                                          "double const* const",
                                          "foo::bar<char> const &",
#ifdef _LIBCPP_VERSION
                                          "foo::bar<foo::bar<char, double>>"};
#else
                                          "foo::bar<foo::bar<char, double> >"};
#endif

public:
    // Returns the name of the type according to the list of names defined above.
    std::string expected_name()
    {
        return names[seqan3::list_traits::find<param_type, as_type_list_t>];
    }
};

// Register test.
TYPED_TEST_SUITE(type_inspection, reflection_types, );

TYPED_TEST(type_inspection, type_name_as_string)
{
    EXPECT_EQ(seqan3::detail::type_name_as_string<TypeParam>, this->expected_name());
}
