// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::band_static.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/detail/strong_type.hpp>
#include <seqan3/std/concepts>

namespace seqan3
{

/*!\brief Type for a lower boundary.
 * \todo Put into core module. This might be useful in other places as well.
 * \ingroup alignment_band
 * \tparam value_t The underlying type of the lower bound; must model seqan3::Arithmetic.
 */
template <seqan3::Arithmetic value_t>
struct lower_bound : detail::strong_type<value_t, lower_bound<value_t>>
{
    //!\brief Inheriting constructors from base class.
    using detail::strong_type<value_t, lower_bound<value_t>>::strong_type;
};

/*!\brief Type for an upper boundary.
 * \todo Put into core module. This might be useful in other places as well.
 * \ingroup alignment_band
 * \tparam value_t The underlying type of the upper bound; must model seqan3::Arithmetic.
 */
template <seqan3::Arithmetic value_t>
struct upper_bound : detail::strong_type<value_t, upper_bound<value_t>>
{
    //!\brief Inheriting constructors from base class.
    using detail::strong_type<value_t, upper_bound<value_t>>::strong_type;
};

/*!\name Deduction guides
 * \brief Deduces template parameter from the argument.
 * \{
 */
/*!\brief Deduces the underlying lower boundary type.
 * \relates seqan3::lower_bound
 * \tparam value_t The underlying type of the lower bound; must model seqan3::Arithmetic.
 */
template <seqan3::Arithmetic value_t>
lower_bound(value_t) -> lower_bound<value_t>;

/*!\brief Deduces the underlying upper boundary type.
 * \relates seqan3::upper_bound
 * \tparam value_t The underlying type of the upper bound; must model seqan3::Arithmetic.
 */
template <seqan3::Arithmetic value_t>
upper_bound(value_t) -> upper_bound<value_t>;
//!\}

/*!\brief Data structure for a static band.
 * \ingroup alignment_band
 */
class static_band
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr static_band()                                noexcept = default; //!< Defaulted
    constexpr static_band(static_band const &)             noexcept = default; //!< Defaulted
    constexpr static_band(static_band &&)                  noexcept = default; //!< Defaulted
    constexpr static_band & operator=(static_band const &) noexcept = default; //!< Defaulted
    constexpr static_band & operator=(static_band &&)      noexcept = default; //!< Defaulted
    ~static_band()                                         noexcept = default; //!< Defaulted

    /*!\brief Construction from seqan3::lower_bound and seqan3::upper_bound.
     * \tparam input_value_t The input type of the lower and upper band boundaries.
     * \param lower The lower boundary of the band; must model std::Integral.
     * \param upper The upper boundary of the band; must model std::Integral.
     *
     * \throws std::invalid_argument if upper < lower.
     *
     * \details
     * The boundaries denote the maximum allowed inbalance of insertions and deletions in the alignment.
     * For a symmetric band, choose lower = -upper. The upper boundary must not be smaller than the lower boundary.
     */
    template <std::Integral input_value_t>
    constexpr static_band(lower_bound<input_value_t> const lower, upper_bound<input_value_t> const upper)
        : lower_bound{lower.get()}, upper_bound{upper.get()}
    {
        if (lower.get() > upper.get())
        {
            throw std::invalid_argument("The upper boundary must not be smaller than the lower boundary.");
        }
    }
    //!}

    //!\brief The data member storing the lower boundary of the band.
    int64_t lower_bound{std::numeric_limits<int64_t>::lowest()};
    //!\brief The data member storing the upper boundary of the band.
    int64_t upper_bound{std::numeric_limits<int64_t>::max()};
};

} // namespace seqan3
