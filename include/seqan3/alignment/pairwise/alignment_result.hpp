// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::alignment_result.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/type_traits/template_inspection.hpp>

namespace seqan3::detail
{

/*!\brief A struct that contains the actual alignment result data.
 * \tparam id_t          The type for the alignment identifier.
 * \tparam score_t       The type for the resulting score.
 * \tparam back_coord_t  The type for the back coordinate, can be omitted.
 * \tparam front_coord_t The type for the front coordinate, can be omitted.
 * \tparam alignment_t   The type for the alignment, can be omitted.
 */
template <typename id_t,
          typename score_t,
          typename back_coord_t = std::nullopt_t *,
          typename front_coord_t = std::nullopt_t *,
          typename alignment_t = std::nullopt_t *>
struct alignment_result_value_type
{
    //! \brief The alignment identifier.
    id_t id{};
    //! \brief The alignment score.
    score_t score{};
    //! \brief The back coordinate of the alignment.
    back_coord_t back_coordinate{};
    //! \brief The front coordinate of the alignment.
    front_coord_t front_coordinate{};
    //! \brief The alignment, i.e. the actual base pair matching.
    alignment_t alignment{};
};

/*!\name Type deduction guides
 * \brief Type deduction for the different combinations of result types.
 * \{
 */
 //! \brief Type deduction for an empty object. It will always fail the compilation, if any field is accessed.
alignment_result_value_type()
    -> alignment_result_value_type<std::nullopt_t *, std::nullopt_t *>;

//! \brief Type deduction for id and score only.
template <typename id_t, typename score_t>
alignment_result_value_type(id_t, score_t)
    -> alignment_result_value_type<id_t, score_t>;

//! \brief Type deduction for id, score and back coordinate.
template <typename id_t, typename score_t, typename back_coord_t>
alignment_result_value_type(id_t, score_t, back_coord_t)
    -> alignment_result_value_type<id_t, score_t, back_coord_t>;

//! \brief Type deduction for id, score, back coordinate and front coordinate.
template <typename id_t, typename score_t, typename back_coord_t, typename front_coord_t>
alignment_result_value_type(id_t, score_t, back_coord_t, front_coord_t)
    -> alignment_result_value_type<id_t, score_t, back_coord_t, front_coord_t>;

//! \brief Type deduction for id, score, back coordinate, front coordinate and alignment.
template <typename id_t, typename score_t, typename back_coord_t, typename front_coord_t, typename alignment_t>
alignment_result_value_type(id_t, score_t, back_coord_t, front_coord_t, alignment_t)
    -> alignment_result_value_type<id_t, score_t, back_coord_t, front_coord_t, alignment_t>;
//!\}

} // namespace seqan3::detail

namespace seqan3
{

/*!\brief Stores the alignment results and gives access to score, alignment and the front and back coordinates.
 * \ingroup pairwise_alignment
 * \tparam alignment_result_traits The type of the traits object.
 *
 * \details
 *
 * Objects of this class are the result of an alignment computation.
 * It always contains an alignment identifier and the resulting score.
 * Optionally – if the user requests – also the begin and end positions within
 * the sequences and the alignment can be calculated. When accessing a field that
 * has not been calculated, an assertion will fail during compilation.
 */
template <typename alignment_result_traits>
//!\cond
    requires detail::is_type_specialisation_of_v<alignment_result_traits, detail::alignment_result_value_type>
//!\endcond
class alignment_result
{
private:
    //! \brief Traits object that contains the actual alignment result data.
    alignment_result_traits data;

    /*!\name Member types
     * \brief Local definition of the types contained in the `data` object.
     * \{
     */
    //! \brief The type for the alignment identifier.
    using id_t          = decltype(data.id);
    //! \brief The type for the resulting score.
    using score_t       = decltype(data.score);
    //! \brief The type for the back coordinate.
    using back_coord_t   = decltype(data.back_coordinate);
    //! \brief The type for the front coordinate.
    using front_coord_t = decltype(data.front_coordinate);
    //! \brief The type for the alignment.
    using alignment_t   = decltype(data.alignment);
    //!\}

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */

    /*!\brief Constructs a seqan3::alignment_result from an `alignment_result_traits` object.
     * \param[in] value The alignment results.
     */
    alignment_result(alignment_result_traits value) : data(value) {};

    alignment_result() = default;                                     //!< Defaulted
    alignment_result(alignment_result const &) = default;             //!< Defaulted
    alignment_result(alignment_result &&) = default;                  //!< Defaulted
    alignment_result & operator=(alignment_result const &) = default; //!< Defaulted
    alignment_result & operator=(alignment_result &&) = default;      //!< Defaulted
    ~alignment_result() = default;                                    //!< Defaulted
    //!\}

    /*!\name Access functions
     * \brief Functions to access elements of the alignment result type.
     * \{
     */

    /*!\brief Returns the alignment identifier.
     * \return The id field.
     * \attention This function will fail the compilation, if the id is not set.
     */
    constexpr id_t id() const noexcept
    {
        static_assert(!std::is_same_v<id_t, std::nullopt_t *>,
                      "Failed to access the identifier.");
        return data.id;
    }

    /*!\brief Returns the alignment score.
     * \return The score field.
     * \attention This function will fail the compilation, if the score is not set.
     */
    constexpr score_t score() const noexcept
    {
        static_assert(!std::is_same_v<score_t, std::nullopt_t *>,
                      "Failed to access the score.");
        return data.score;
    }

    /*!\brief Returns the back coordinate of the alignment.
     * \return A pair of positions in the respective sequences, where the calculated alignment ends (inclusive).
     * \attention This function will fail the compilation, if the back coordinate was not requested in the alignment
     * configuration.
     */
    constexpr back_coord_t const & back_coordinate() const noexcept
    {
        static_assert(!std::is_same_v<back_coord_t, std::nullopt_t *>,
                      "Trying to access the back coordinate, although it was not requested in the alignment "
                      "configuration.");
        return data.back_coordinate;
    }

    /*!\brief Returns the front coordinate of the alignment.
     * \return  A pair of positions in the respective sequences, where the calculated alignment starts.
     * \details Guaranteed to be smaller than or equal to `back_coordinate()`.
     * \attention This function will fail the compilation, if the front coordinate was not requested in the alignment
     * configuration.
     */
    constexpr front_coord_t const & front_coordinate() const noexcept
    {
        static_assert(!std::is_same_v<front_coord_t, std::nullopt_t *>,
                      "Trying to access the front coordinate, although it was not requested in the alignment "
                      "configuration.");
        return data.front_coordinate;
    }

    /*!\brief Returns the actual alignment, i.e. the base pair matching.
     * \return At least two aligned sequences, which represent the alignment.
     * \attention This function with fail the compilation, if the alignment was not requested in the alignment
     * configuration.
     */
    constexpr alignment_t const & alignment() const noexcept
    {
        static_assert(!std::is_same_v<alignment_t, std::nullopt_t *>,
                      "Trying to access the alignment, although it was not requested in the alignment configuration.");
        return data.alignment;
    }
    //!\}
};

} // namespace seqan3
