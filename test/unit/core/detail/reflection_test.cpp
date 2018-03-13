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
using reflection_types = ::testing::Types<char, char16_t, char32_t,
                                          short int, unsigned int, double, const char *,
                                          foo::bar<char>, foo::bar<foo::bar<char, double>>>;

template <typename param_type>
class reflection : public ::testing::Test
{
    // The corresponding list of names that should be generated. Must have the same order as `reflection_types`.
    inline static const std::vector names{"char", "char16_t", "char32_t", "short int", "unsigned int", "double",
                                          "const char*", "foo::bar<char>", "foo::bar<foo::bar<char, double> >"};

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

TYPED_TEST(reflection, name)
{
    EXPECT_EQ(detail::get_display_name_v<TypeParam>.string(), this->expected_name());
}
