// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \brief Provides the data structures and precomputed instances for (optimum) search schemes.
 */

#pragma once

#include <array>
#include <vector>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

//!\brief Object storing information for a search (of a search scheme).
//!\attention Number of blocks have to be known at compile time.
//!\ingroup search
template <uint8_t nbr_blocks>
struct search
{
    //!\brief Type for storing the length of blocks
    using blocks_length_type = std::array<size_t, nbr_blocks>;

    //!\brief Order of blocks
    std::array<uint8_t, nbr_blocks> pi;
    //!\brief Lower error bound for each block (accumulated values)
    std::array<uint8_t, nbr_blocks> l;
    //!\brief Upper error bound for each block (accumulated values)
    std::array<uint8_t, nbr_blocks> u;

    //!\brief Returns the number of blocks
    constexpr uint8_t blocks() const noexcept
    {
        return nbr_blocks;
    }
};

//!\brief Object storing information for a search (of a search scheme).
//!\attention Number of blocks do not have to be known at compile time.
//!\ingroup search
struct search_dyn
{
    //!\brief Type for storing the length of blocks
    using blocks_length_type = std::vector<size_t>;

    //!\brief Order of blocks
    std::vector<uint8_t> pi;
    //!\brief Lower error bound for each block (accumulated values)
    std::vector<uint8_t> l;
    //!\brief Upper error bound for each block (accumulated values)
    std::vector<uint8_t> u;

    //!\brief Returns the number of blocks
    uint8_t blocks() const noexcept
    {
        return pi.size();
    }
};

//!\brief Type for storing search schemes. Number of blocks have to be known at compile time.
//!\ingroup search
template <uint8_t nbr_searches, uint8_t nbr_blocks>
using search_scheme_type = std::array<search<nbr_blocks>, nbr_searches>;

//!\brief Type for storing search schemes. Number of blocks do not have to be known at compile time.
//!\ingroup search
using search_scheme_dyn_type = std::vector<search_dyn>;

/*!\brief Search scheme that is optimal in the running time for the specified lower and upper error bound.
 * \ingroup search
 * \tparam min_error Lower bound of errors.
 * \tparam max_error Upper bound of errors.
 * \details Please note that the searches within each search scheme are sorted by their asymptotical run time
 *          (i.e. upper error bound string), s.t. easy to compute searches come first. This improves the run time of
 *          algorithms that abort after the first hit (e.g. search mode: best). Even though it is not guaranteed, this
 *          seems to be a good greedy approach.
 */
template <uint8_t min_error, uint8_t max_error>
inline constexpr int optimum_search_scheme{0};

//!\cond

template <>
inline constexpr search_scheme_type<1, 1> optimum_search_scheme<0, 0>{{{{1}, {0}, {0}}}};

template <>
inline constexpr search_scheme_type<2, 2> optimum_search_scheme<0, 1>{
    {{{1, 2}, {0, 0}, {0, 1}}, {{2, 1}, {0, 1}, {0, 1}}}};

template <>
inline constexpr search_scheme_type<2, 2> optimum_search_scheme<1, 1>{
    {{{1, 2}, {0, 1}, {0, 1}}, {{2, 1}, {0, 1}, {0, 1}}}};

template <>
inline constexpr search_scheme_type<3, 4> optimum_search_scheme<0, 2>{{{{1, 2, 3, 4}, {0, 0, 1, 1}, {0, 0, 2, 2}},
                                                                       {{3, 2, 1, 4}, {0, 0, 0, 0}, {0, 1, 1, 2}},
                                                                       {{4, 3, 2, 1}, {0, 0, 0, 2}, {0, 1, 2, 2}}}};

template <>
inline constexpr search_scheme_type<3, 4> optimum_search_scheme<1, 2>{{{{1, 2, 3, 4}, {0, 0, 0, 1}, {0, 0, 2, 2}},
                                                                       {{3, 2, 1, 4}, {0, 0, 1, 1}, {0, 1, 1, 2}},
                                                                       {{4, 3, 2, 1}, {0, 0, 0, 2}, {0, 1, 2, 2}}}};

template <>
inline constexpr search_scheme_type<3, 4> optimum_search_scheme<2, 2>{{{{4, 3, 2, 1}, {0, 0, 1, 2}, {0, 0, 2, 2}},
                                                                       {{2, 3, 4, 1}, {0, 0, 0, 2}, {0, 1, 1, 2}},
                                                                       {{1, 2, 3, 4}, {0, 0, 0, 2}, {0, 1, 2, 2}}}};

// TODO: benchmark whether the first search is really the fastest one (see \details of optimum_search_scheme)
template <>
inline constexpr search_scheme_type<4, 5> optimum_search_scheme<0, 3>{
    {{{5, 4, 3, 2, 1}, {0, 0, 0, 0, 0}, {0, 0, 3, 3, 3}},
     {{3, 4, 5, 2, 1}, {0, 0, 1, 1, 1}, {0, 1, 1, 2, 3}},
     {{2, 3, 4, 5, 1}, {0, 0, 0, 2, 2}, {0, 1, 2, 2, 3}},
     {{1, 2, 3, 4, 5}, {0, 0, 0, 0, 3}, {0, 2, 2, 3, 3}}}};

template <>
inline constexpr search_scheme_type<4, 5> optimum_search_scheme<1, 3>{
    {{{5, 4, 3, 2, 1}, {0, 0, 0, 0, 1}, {0, 0, 3, 3, 3}},
     {{3, 4, 5, 2, 1}, {0, 0, 1, 1, 1}, {0, 1, 1, 2, 3}},
     {{2, 3, 4, 5, 1}, {0, 0, 0, 2, 2}, {0, 1, 2, 2, 3}},
     {{1, 2, 3, 4, 5}, {0, 0, 0, 0, 3}, {0, 2, 2, 3, 3}}}};

template <>
inline constexpr search_scheme_type<4, 5> optimum_search_scheme<2, 3>{
    {{{5, 4, 3, 2, 1}, {0, 0, 0, 0, 2}, {0, 0, 3, 3, 3}},
     {{3, 4, 5, 2, 1}, {0, 0, 1, 1, 2}, {0, 1, 1, 2, 3}},
     {{2, 3, 4, 5, 1}, {0, 0, 0, 2, 2}, {0, 1, 2, 2, 3}},
     {{1, 2, 3, 4, 5}, {0, 0, 0, 0, 3}, {0, 2, 2, 3, 3}}}};

template <>
inline constexpr search_scheme_type<4, 5> optimum_search_scheme<3, 3>{
    {{{5, 4, 3, 2, 1}, {0, 0, 0, 0, 3}, {0, 0, 3, 3, 3}},
     {{3, 4, 5, 2, 1}, {0, 0, 1, 1, 3}, {0, 1, 1, 2, 3}},
     {{2, 3, 4, 5, 1}, {0, 0, 0, 2, 3}, {0, 1, 2, 2, 3}},
     {{1, 2, 3, 4, 5}, {0, 0, 0, 0, 3}, {0, 2, 2, 3, 3}}}};

// TODO: add the following missing optimum search schemes (computation has not finished yet)
// optimum_search_scheme<i, 4>, 0 < i <= 4

//!\endcond

} // namespace seqan3::detail
