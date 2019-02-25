// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::align_result.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/metafunction/template_inspection.hpp>

namespace seqan3::detail
{

/*!\brief A struct that contains the actual alignment result data.
 * \tparam id_t          The type for the alignment identifier.
 * \tparam score_t       The type for the resulting score.
 * \tparam end_coord_t   The type for the end coordinate, can be omitted.
 * \tparam begin_coord_t The type for the begin coordinate, can be omitted.
 * \tparam alignment_t   The type for the alignment, can be omitted.
 */
template <typename id_t,
          typename score_t,
          typename end_coord_t = std::nullopt_t *,
          typename begin_coord_t = std::nullopt_t *,
          typename alignment_t = std::nullopt_t *>
struct align_result_value_type
{
    //! \brief The alignment identifier.
    id_t id{};
    //! \brief The alignment score.
    score_t score{};
    //! \brief The end coordinate of the alignment.
    end_coord_t end_coordinate{};
    //! \brief The begin coordinate of the alignment.
    begin_coord_t begin_coordinate{};
    //! \brief The alignment, i.e. the actual base pair matching.
    alignment_t alignment{};
};

/*!\name Type deduction guides
 * \brief Type deduction for the different combinations of result types.
 * \{
 */
 //! \brief Type deduction for an empty object. It will always fail the compilation, if any field is accessed.
align_result_value_type()
    -> align_result_value_type<std::nullopt_t *, std::nullopt_t *>;

//! \brief Type deduction for id and score only.
template <typename id_t, typename score_t>
align_result_value_type(id_t, score_t)
    -> align_result_value_type<id_t, score_t>;

//! \brief Type deduction for id, score and end coordinate.
template <typename id_t, typename score_t, typename end_coord_t>
align_result_value_type(id_t, score_t, end_coord_t)
    -> align_result_value_type<id_t, score_t, end_coord_t>;

//! \brief Type deduction for id, score, end coordinate and begin coordinate.
template <typename id_t, typename score_t, typename end_coord_t, typename begin_coord_t>
align_result_value_type(id_t, score_t, end_coord_t, begin_coord_t)
    -> align_result_value_type<id_t, score_t, end_coord_t, begin_coord_t>;

//! \brief Type deduction for id, score, end coordinate, begin coordinate and alignment.
template <typename id_t, typename score_t, typename end_coord_t, typename begin_coord_t, typename alignment_t>
align_result_value_type(id_t, score_t, end_coord_t, begin_coord_t, alignment_t)
    -> align_result_value_type<id_t, score_t, end_coord_t, begin_coord_t, alignment_t>;
//!\}

} // namespace seqan3::detail

namespace seqan3
{

/*!\brief Stores the alignment results and gives access to score, alignment and the begin and end coordinates.
 * \ingroup pairwise
 * \tparam align_result_traits The type of the traits object.
 *
 * \details
 *
 * Objects of this class are the result of an alignment computation.
 * It always contains an alignment identifier and the resulting score.
 * Optionally – if the user requests – also the begin and end positions within
 * the sequences and the alignment can be calculated. When accessing a field that
 * has not been calculated, an assertion will fail during compilation.
 */
template <typename align_result_traits>
//!\cond
    requires detail::is_type_specialisation_of_v<align_result_traits, detail::align_result_value_type>
//!\endcond
class align_result
{
private:
    //! \brief Traits object that contains the actual alignment result data.
    align_result_traits data;

    /*!\name Member types
     * \brief Local definition of the types contained in the `data` object.
     * \{
     */
    //! \brief The type for the alignment identifier.
    using id_t          = decltype(data.id);
    //! \brief The type for the resulting score.
    using score_t       = decltype(data.score);
    //! \brief The type for the end coordinate.
    using end_coord_t   = decltype(data.end_coordinate);
    //! \brief The type for the begin coordinate.
    using begin_coord_t = decltype(data.begin_coordinate);
    //! \brief The type for the alignment.
    using alignment_t   = decltype(data.alignment);
    //!\}

public:
    /*!\name Constructors, destructor and assignment
     * \{
     *
     * \brief Constructor to pass the alignment result traits.
     * \param[in] value The alignment results.
     */
    align_result(align_result_traits value) : data(value) {};

    //! \brief Default constructor.
    align_result() = default;
    //! \brief Default copy constructor.
    align_result(align_result const &) = default;
    //! \brief Default move constructor.
    align_result(align_result &&) = default;
    //! \brief Default copy assignment.
    align_result & operator=(align_result const &) = default;
    //! \brief Default move assignment.
    align_result & operator=(align_result &&) = default;
    //! \brief Default destructor.
    ~align_result() = default;
    //!\}

    /*!\name Access functions
     * \brief Functions to access elements of the alignment result type.
     * \{
     */

    /*!\brief Returns the alignment identifier.
     * \return The id field.
     * \attention This function with fail the compilation, if the id is not set.
     */
    constexpr id_t get_id() const noexcept
    {
        static_assert(!std::is_same_v<id_t, std::nullopt_t *>,
                      "Failed to access the identifier.");
        return data.id;
    }

    /*!\brief Returns the alignment score.
     * \return The score field.
     * \attention This function with fail the compilation, if the score is not set.
     */
    constexpr score_t get_score() const noexcept
    {
        static_assert(!std::is_same_v<score_t, std::nullopt_t *>,
                      "Failed to access the score.");
        return data.score;
    }

    /*!\brief Returns the end coordinate of the alignment.
     * \return A pair of positions in the respective sequences, where the calculated alignment ends.
     * \attention This function with fail the compilation, if the end coordinate was not requested in the alignment
     * configuration.
     */
    constexpr end_coord_t const & get_end_coordinate() const noexcept
    {
        static_assert(!std::is_same_v<end_coord_t, std::nullopt_t *>,
                      "Trying to access the end coordinate, although it was not requested in the alignment "
                      "configuration.");
        return data.end_coordinate;
    }

    /*!\brief Returns the begin coordinate of the alignment.
     * \return  A pair of positions in the respective sequences, where the calculated alignment starts.
     * \details Guaranteed to be smaller than or equal to `get_end_coordinate()`.
     * \attention This function with fail the compilation, if the begin coordinate was not requested in the alignment
     * configuration.
     */
    constexpr begin_coord_t const & get_begin_coordinate() const noexcept
    {
        static_assert(!std::is_same_v<begin_coord_t, std::nullopt_t *>,
                      "Trying to access the begin coordinate, although it was not requested in the alignment "
                      "configuration.");
        return data.begin_coordinate;
    }

    /*!\brief Returns the actual alignment, i.e. the base pair matching.
     * \return At least two gapped sequences, which represent the alignment.
     * \attention This function with fail the compilation, if the alignment was not requested in the alignment
     * configuration.
     */
    constexpr alignment_t const & get_alignment() const noexcept
    {
        static_assert(!std::is_same_v<alignment_t, std::nullopt_t *>,
                      "Trying to access the alignment, although it was not requested in the alignment configuration.");
        return data.alignment;
    }
    //!\}
};

} // namespace seqan3
