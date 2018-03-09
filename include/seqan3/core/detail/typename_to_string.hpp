// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
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

#pragma once

#include <type_traits>
#include <typeinfo>

#if __has_include(<cxxabi.h>)
#include <memory>
#include <cxxabi.h>
#endif

#include <seqan3/core/platform.hpp>

/*!\file
 * \brief Utility function to extract the human readable name of a type using typeinfo
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \ingroup core
 */

namespace seqan3::detail
{

//!\cond DEV
/*!\brief Given an instance of \a entity_t returns the human readable name of this type.
 * \tparam entity_t The type to query it's name for.
 * \param entity An instance of \a entity_t.
 *
 * This function uses typeid to query the name of the type.
 * On some platforms the returned name is mangled.
 * In this case the function tries to demangle the name.
 * In case of an failure it returns the mangled name.
 *
 * \returns std::string with the name of the type.
 *
 * \par Exception
 * Might throw either std::bad_typeid or std::bad_alloc in creating the std::string.
 *
 * \par Example
 * \code{.cpp}
 * #include <iostream>
 * #include <seqan3/core/detail/typename_to_string.hpp>
 *
 * template <typename type = void>
 * struct foo
 * {};
 *
 * int main()
 * {
 *     std::cout << "The type's name is: " << seqan3::detail::typename_to_string(foo<>{}) << std::endl;
 *     return 0;
 * }
 * \endcode
 * The above example will produce the following output:
 * \code{.sh}
 * The type's name is: foo<void>
 * \endcode
 */
template <typename entity_t>
inline std::string
typename_to_string(entity_t const & entity)
{
#if __has_include(<cxxabi.h>)  // We should use demangle.
    int status;

    auto deleter = [] (char * rsrc) { free(rsrc); };
    std::unique_ptr<char, decltype(deleter)> name_rsrc{
        abi::__cxa_demangle(typeid(entity).name(), 0, 0, &status),
        deleter};

    if (status == 0)
        return {name_rsrc.get()};
#endif // __has_include(<cxxabi.h>)
    return {typeid(entity).name()};
}
//!\endcond

} // namespace seqan3::detail
