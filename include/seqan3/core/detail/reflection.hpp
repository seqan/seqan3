// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \brief Static reflection for arbitrary types.
 */

#pragma once

#include <utility>

#include <seqan3/range/container/small_string.hpp>

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
 * \include test/snippet/core/detail/reflection.cpp
 */
template <typename type>
struct get_display_name
{
private:
    //!\brief Helper function to get the display name.
    static constexpr auto get_display_name_fn()
    {
        // Use a helper function to extract the size of the requested type.
        constexpr auto name_length = get_display_name_size_v<type>;

        // Extract the type again.
        auto name_ptr = __PRETTY_FUNCTION__;

        // Move pointer to first letter of actual type name.
        while (*name_ptr++ != '=')
        {}

        for (; *name_ptr == ' '; ++name_ptr)
        {}

        return small_string<name_length>{name_ptr, name_ptr + name_length};
    };
public:
    //!\brief Defines the display name of `type`.
    static constexpr small_string value = get_display_name_fn();
};

//!\brief Shortcut for accessing the display name of the given `type`.
//!\see get_display_name
//!\ingroup core
template <typename type>
constexpr small_string get_display_name_v = get_display_name<type>::value;
} // namespace seqan3::detail
