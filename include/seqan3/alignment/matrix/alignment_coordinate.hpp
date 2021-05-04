// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Provides seqan3::alignment_coordinate.
 */

#pragma once

#include <seqan3/alignment/matrix/detail/advanceable_alignment_coordinate.hpp>
#include <seqan3/alignment/matrix/detail/matrix_coordinate.hpp>
#include <seqan3/core/detail/debug_stream_tuple.hpp>
#include <seqan3/core/detail/debug_stream_type.hpp>

namespace seqan3
{

#ifdef SEQAN3_DEPRECATED_310
/*!\brief Represents the begin/end of the pairwise alignment in the respective sequences.
 * \ingroup alignment_matrix
 *
 * \if DEV
 * \details
 * This class only gives access to the respective positions of the sequences and is meant for
 * the user interface. The additional complexity of an advanceable coordinate using the
 * seqan3::detail::advanceable_alignment_coordinate is only necessary for the implementation of the pairwise
 * alignment algorithm. Within in the algorithm the coordinate is used in combination with a seqan3::views::iota to
 * keep track of the current position within the alignment matrix. For the user, however, this interface adds no
 * benefit as they are only interested in the front/back coordinates for the respective alignment.
 * \endif
 *
 * \deprecated This class will be removed in SeqAn 3.1.
 */
class SEQAN3_DEPRECATED_310 alignment_coordinate
//!\cond DEV
    : public detail::advanceable_alignment_coordinate<detail::advanceable_alignment_coordinate_state::none>
//!\endcond
{
    //!\cond DEV
    //!\brief The type of the base class.
    using base_t = detail::advanceable_alignment_coordinate<detail::advanceable_alignment_coordinate_state::none>;
    //!\endcond

public:

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr alignment_coordinate() = default;                                         //!< Defaulted
    constexpr alignment_coordinate(alignment_coordinate const &) = default;             //!< Defaulted
    constexpr alignment_coordinate(alignment_coordinate &&) = default;                  //!< Defaulted
    constexpr alignment_coordinate & operator=(alignment_coordinate const &) = default; //!< Defaulted
    constexpr alignment_coordinate & operator=(alignment_coordinate &&) = default;      //!< Defaulted
    ~alignment_coordinate() = default;                                                  //!< Defaulted

    //!\cond DEV
    //!\brief Inherit the constructor from the base class.
    using base_t::base_t;

    //!\brief Constructs from the seqan3::detail::advanceable_alignment_coordinate base class.
    constexpr alignment_coordinate(base_t const & base) : base_t{base}
    {}

    //!\brief Constructs from the seqan3::detail::advanceable_alignment_coordinate base class.
    constexpr alignment_coordinate(base_t && base) : base_t{std::move(base)}
    {}
    //!\endcond
    //!\}

    using base_t::first;
    using base_t::second;

    //!\brief The begin/end position of the alignment in the first sequence.
    SEQAN3_DOXYGEN_ONLY(size_t first;)
    //!\brief The begin/end position of the alignment in the second sequence.
    SEQAN3_DOXYGEN_ONLY(size_t second;)

    //!\privatesection
    //!\brief Implicit conversion to seqan3::detail::matrix_coordinate.
    constexpr operator detail::matrix_coordinate() const
    {
        return detail::matrix_coordinate{detail::row_index_type{second}, detail::column_index_type{first}};
    }
};

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
//!\cond
template <typename char_t, typename coordinate_type>
   requires std::same_as<std::remove_cvref_t<coordinate_type>, alignment_coordinate>
inline debug_stream_type<char_t> & operator<<(debug_stream_type<char_t> & s, coordinate_type && c)
{
   s << std::tie(c.first, c.second);
   return s;
}
//!\endcond
#pragma GCC diagnostic pop
#endif // SEQAN3_DEPRECATED_310

} // namespace seqan3
