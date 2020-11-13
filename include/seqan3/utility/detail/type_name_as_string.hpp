// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides traits to inspect some information of a type, for example its name.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#if defined(__GNUC__) || defined(__clang__)
#include <cxxabi.h>
#endif // defined(__GNUC__) || defined(__clang__)

#include <functional>
#include <memory>
#include <string>
#include <typeinfo>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

/*!\brief Defines the human-readable name of the given type using the
          [typeid](https://en.cppreference.com/w/cpp/language/typeid) operator.
 * \ingroup core
 *
 * \tparam type The type to get the human-readable name for.
 *
 * \details
 *
 * On gcc and clang std::type_info only returns a mangled name.
 * The mangled name can be converted to human-readable form using implementation-specific API such as
 * abi::__cxa_demangle. In other implementations the name returned is already human-readable.
 *
 * \note The returned name is implementation defined and might change between different tool chains.
 */
template <typename type>
inline std::string const type_name_as_string = [] ()
{
    std::string demangled_name{};
#if defined(__GNUC__) || defined(__clang__) // clang and gcc only return a mangled name.
    using safe_ptr_t = std::unique_ptr<char, std::function<void(char *)>>;

    int status{};
    safe_ptr_t demangled_name_ptr{abi::__cxa_demangle(typeid(type).name(), 0, 0, &status),
                                  [] (char * name_ptr) { free(name_ptr); }};
    demangled_name = std::string{std::addressof(*demangled_name_ptr)};
#else // e.g. MSVC
    demangled_name = typeid(type).name();
#endif // defined(__GNUC__) || defined(__clang__)

    if constexpr (std::is_const_v<std::remove_reference_t<type>>)
        demangled_name += " const";
    if constexpr (std::is_lvalue_reference_v<type>)
        demangled_name += " &";
    if constexpr (std::is_rvalue_reference_v<type>)
        demangled_name += " &&";

    return demangled_name;
}();

}  // namespace seqan3::detail
