// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
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
// ============================================================================

/*!\file
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \brief Static reflection for arbitrary types.
 */

#pragma once

#include <utility>

#include <seqan3/range/container/constexpr_string.hpp>

namespace seqan3::detail
{

/*!\brief Extracts the size of the display name of the given template.
 * \tparam type The template parameter to get the size of the display name for.
 * \ingroup core
 * \see get_display_name_size_v
 *
 * This functions uses gcc's [__PRETTY_FUNCTION__](https://gcc.gnu.org/onlinedocs/gcc/Function-Names.html)
 * to get the size of the display name.
 * This code adapts the snippet from the following
 * [discussion](https://stackoverflow.com/questions/35941045/can-i-obtain-c-type-names-in-a-constexpr-way/35943472#35943472).
 */
template <typename type>
struct get_display_name_size
{
private:
    //!\brief Helper function to extract the size of the display name.
    static constexpr auto get_size()
    {
        // __PRETTY_FUNCTION__ exposes the signature of the class including the name of the template instance as
        // a static const char[]. The following code, extracts the part that displays the template name.
        auto name_ptr = __PRETTY_FUNCTION__;

        // Move pointer to first letter of actual type name.
        while (*name_ptr++ != '=')
        {}

        for (; *name_ptr == ' '; ++name_ptr)
        {}

        // Find the end of the actual type name.
        char const * end_name_ptr = name_ptr;
        int count = 1;

        for (; ; ++end_name_ptr)
        {
            switch (*end_name_ptr)
            {
                case '[':
                    ++count;
                    break;
                case ']':
                    --count;
                    if (!count)
                        return size_t(end_name_ptr - name_ptr);
                    break;
            }
        }
        return size_t(0);
    };

public:
    //!\brief Value constant storing the size of the display name of `type`.
    static constexpr size_t value = get_size();
};

//!\brief Shortcut for accessing the size of the display name of `type`.
//!\see get_display_name_size
//!\ingroup core
template <typename type>
constexpr size_t get_display_name_size_v = get_display_name_size<type>::value;

/*!\brief Metafunction to get the display name of a given type.
 * \tparam type The type to get the display name for.
 * \ingroup core
 * \see get_display_name_v
 *
 * Uses gcc's [__PRETTY_FUNCTION__](https://gcc.gnu.org/onlinedocs/gcc/Function-Names.html) to get the
 * function signature including the template name.
 *
 * ### Example
 *
 * The following snippet demonstrates the usage:
 * ```cpp
 * #include <iostream>
 *
 * #include <seqan3/core/detail/reflection.hpp>
 *
 * namespace foo
 * {
 *
 * template <typename ...type>
 * struct bar
 * {};
 *
 * } // namespace foo
 *
 * int main()
 * {
 *     std::cout << detail::get_display_name_v<foo::bar<char, double>>.string() << std::endl; // prints: foo::bar<char, double> >
 * }
 * ```
 */
template <typename type>
struct get_display_name
{
private:
    //!\brief Helper function to get the display name.
    static constexpr auto get_display_name_fn()
    {
        // Use a helper function to extract the size of the requested type.
        constexpr_string<get_display_name_size_v<type>> tmp{};

        // Extract the type again.
        auto name_ptr = __PRETTY_FUNCTION__;

        // Move pointer to first letter of actual type name.
        while (*name_ptr++ != '=')
        {}

        for (; *name_ptr == ' '; ++name_ptr)
        {}

        for (unsigned i = 0; i < tmp.size(); ++i, ++name_ptr)
        {
            tmp[i] = *name_ptr;
        }
        return tmp;
    };
public:
    //!\brief Defines the display name of `type`.
    static constexpr constexpr_string value = get_display_name_fn();
};

//!\brief Shortcut for accessing the display name of the given `type`.
//!\see get_display_name
//!\ingroup core
template <typename type>
constexpr constexpr_string get_display_name_v = get_display_name<type>::value;
} // namespace seqan3::detail
