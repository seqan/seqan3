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
 * \brief Resolves dependencies on the concepts.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <functional>
#include <type_traits>

#include <meta/meta.hpp>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{
/*!\interface seqan3::detail::config_concept <>
 * \ingroup core_algorithm
 * \brief Concept for an algorithm configuration.
 */
/*!\name Requirements for seqan3::detail::config_concept
 * \relates seqan3::detail::config_concept
 * \brief You can expect these member types and the free function on all types that satisfy seqan3::detail::config_concept.
 * \{
 */

/*!\fn      seqan3::detail::config_concept::data();
 * \brief   Gives access to the stored configuration state.
 * \returns The stored configuration state.
 *
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */

//!\}
//!\cond
template <typename config_t>
concept bool config_concept = requires (config_t & cfg, config_t const & cfg_c)
{
    { cfg.data() };
    { cfg_c.data() };
    { std::move(cfg).data() };
    { std::move(cfg_c).data() };

    requires std::is_lvalue_reference_v<decltype(cfg.data())>;
    requires std::is_lvalue_reference_v<decltype(cfg_c.data())> &&
             std::is_const_v<std::remove_reference_t<decltype(cfg_c.data())>>;
    requires std::is_rvalue_reference_v<decltype(std::move(cfg).data())>;
    requires std::is_rvalue_reference_v<decltype(std::move(cfg_c).data())> &&
             std::is_const_v<std::remove_reference_t<decltype(std::move(cfg_c).data())>>;
};
//!\endcond

//!\cond
template <detail::config_concept ... configs_t>
    requires sizeof...(configs_t) >= 1
class configurator;
//!\endcond
} // namespace seqan3

namespace seqan3
{
/*!\name Tuple-like get interface
 * \relates seqan3::detail::configurator
 * \{
 */

//!\brief Provides tuple-like get interface.
template <size_t elem_no, typename ... configs_t>
constexpr auto & get(seqan3::detail::configurator<configs_t...> & cfg) noexcept
{
    return cfg.template get<elem_no>();
}

//!\brief Provides tuple-like get interface.
template <size_t elem_no, typename ... configs_t>
constexpr auto const & get(seqan3::detail::configurator<configs_t...> const & cfg) noexcept
{
    return cfg.template get<elem_no>();
}

//!\brief Provides tuple-like get interface.
template <size_t elem_no, typename ... configs_t>
constexpr auto && get(seqan3::detail::configurator<configs_t...> && cfg) noexcept
{
    return std::move(cfg).template get<elem_no>();
}

//!\brief Provides tuple-like get interface.
template <size_t elem_no, typename ... configs_t>
constexpr auto const && get(seqan3::detail::configurator<configs_t...> const && cfg) noexcept
{
    return std::move(cfg).template get<elem_no>();
}

//!\brief Provides tuple-like get interface.
template <typename target_t, typename ... configs_t>
constexpr auto & get(seqan3::detail::configurator<configs_t...> & cfg) noexcept
{
    return cfg.template get<target_t>();
}

//!\brief Provides tuple-like get interface.
template <typename target_t, typename ... configs_t>
constexpr auto const & get(seqan3::detail::configurator<configs_t...> const & cfg) noexcept
{
    return cfg.template get<target_t>();
}

//!\brief Provides tuple-like get interface.
template <typename target_t, typename ... configs_t>
constexpr auto && get(seqan3::detail::configurator<configs_t...> && cfg) noexcept
{
    return std::move(cfg).template get<target_t>();
}

//!\brief Provides tuple-like get interface.
template <typename target_t, typename ... configs_t>
constexpr auto const && get(seqan3::detail::configurator<configs_t...> const && cfg) noexcept
{
    return std::move(cfg).template get<target_t>();
}
//!\}
} // namespace seqan3

namespace std
{

/*!\brief Value metafunction specialisation for seqan3::detail::configurator
 * \ingroup core_algorithm
 * \relates seqan3::detail::configurator
 * \returns The number of configurations contained in seqan3::detail::configurator.
 */
template <typename ... configs_t>
struct tuple_size<seqan3::detail::configurator<configs_t...>>
{
    //!\brief The value member.
    static constexpr size_t value = sizeof...(configs_t);
};

/*!\brief Value metafunction specialisation for seqan3::detail::configurator.
 * \ingroup core_algorithm
 * \relates seqan3::detail::configurator
 * \returns The the type of configuration for the requested position.
 */
template <size_t elem_no, typename ... configs_t>
struct tuple_element<elem_no, seqan3::detail::configurator<configs_t...>>
{
    //!\brief The member type.
    using type = meta::at_c<typename seqan3::detail::configurator<configs_t...>::type_list_type, elem_no>;;
};

//!\cond
//!\brief Provides std::get interface.
template <size_t elem_no, typename ... configs_t>
constexpr auto & get(seqan3::detail::configurator<configs_t...> & cfg) noexcept
{
    return seqan3::get<elem_no>(cfg);
}

//!\brief Provides std::get interface.
template <size_t elem_no, typename ... configs_t>
constexpr auto const & get(seqan3::detail::configurator<configs_t...> const & cfg) noexcept
{
    return  seqan3::get<elem_no>(cfg);
}

//!\brief Provides std::get interface.
template <size_t elem_no, typename ... configs_t>
constexpr auto && get(seqan3::detail::configurator<configs_t...> && cfg) noexcept
{
    return seqan3::get<elem_no>(std::move(cfg));
}

//!\brief Provides std::get interface.
template <size_t elem_no, typename ... configs_t>
constexpr auto const && get(seqan3::detail::configurator<configs_t...> const && cfg) noexcept
{
    return seqan3::get<elem_no>(std::move(cfg));
}

//!\brief Provides std::get interface.
template <typename target_t, typename ... configs_t>
constexpr auto & get(seqan3::detail::configurator<configs_t...> & cfg) noexcept
{
    return seqan3::get<target_t>(cfg);
}

//!\brief Provides std::get interface.
template <typename target_t, typename ... configs_t>
constexpr auto const & get(seqan3::detail::configurator<configs_t...> const & cfg) noexcept
{
    return  seqan3::get<target_t>(cfg);
}

//!\brief Provides std::get interface.
template <typename target_t, typename ... configs_t>
constexpr auto && get(seqan3::detail::configurator<configs_t...> && cfg) noexcept
{
    return seqan3::get<target_t>(std::move(cfg));
}

//!\brief Provides std::get interface.
template <typename target_t, typename ... configs_t>
constexpr auto const && get(seqan3::detail::configurator<configs_t...> const && cfg) noexcept
{
    return seqan3::get<target_t>(std::move(cfg));
}
//!\endcond
} // namespace std
