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
 * \brief Provides seqan3::align_result.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 *
 * \details
 * Objects of this class are the result of an alignment computation.
 * It always contains and alignment identifier and the resulting score.
 * Optionally – if the user requests – also the begin and end positions within
 * the sequences and the trace can be calculated. When accessing a field that
 * has not been calculated an std::bad_optional_access exception is thrown.
 */

#pragma once

#include <optional>

namespace seqan3
{

/*!\brief Stores the alignment results and gives access to score, traceback and the begin and end coordinates.
 * \ingroup pairwise
 * \tparam id_t         The type for the alignment identifier.
 * \tparam score_t      The type for the resulting score.
 * \tparam coordinate_t The type for both the begin and end coordinates, can be omitted.
 * \tparam trace_t      The type for the alignment trace, can be omitted.
 */
template <typename id_t, typename score_t, typename coordinate_t = bool, typename trace_t = bool>
class align_result
{
private:
    //! \brief The alignment identifier.
    id_t id;
    //! \brief The alignment score.
    score_t score;
    //! \brief The end coordinate of the alignment.
    std::optional<coordinate_t> end_coordinate;
    //! \brief The begin coordinate of the alignment.
    std::optional<coordinate_t> begin_coordinate;
    //! \brief The alignment trace, i.e. the actual base pair matching.
    std::optional<trace_t> trace;

public:
    /*!\brief Constructor with all available fields.
     * \param id_               The alignment identifier.
     * \param score_            The alignment score.
     * \param end_coordinate_   The end coordinate of the alignment.
     * \param begin_coordinate_ The begin coordinate of the alignment.
     * \param trace_            The alignment trace, i.e. the actual base pair matching.
     */
    align_result(id_t const id_,
                 score_t const score_,
                 coordinate_t const end_coordinate_,
                 coordinate_t const begin_coordinate_,
                 trace_t const & trace_)
        : id{id_}, score{score_}, end_coordinate{end_coordinate_}, begin_coordinate{begin_coordinate_}, trace{trace_}
    {}

    /*!\brief Constructor without trace.
     * \param id_               The alignment identifier.
     * \param score_            The alignment score.
     * \param end_coordinate_   The end coordinate of the alignment.
     * \param begin_coordinate_ The begin coordinate of the alignment.
     */
    align_result(id_t const id_,
                 score_t const score_,
                 coordinate_t const end_coordinate_,
                 coordinate_t const begin_coordinate_)
        : id{id_}, score{score_}, end_coordinate{end_coordinate_}, begin_coordinate{begin_coordinate_}, trace{}
    {}

    /*!\brief Constructor without trace and begin coordinate information.
     * \param id_               The alignment identifier.
     * \param score_            The alignment score.
     * \param end_coordinate_   The end coordinate of the alignment.
     */
    align_result(id_t const id_,
                 score_t const score_,
                 coordinate_t const end_coordinate_)
        : id{id_}, score{score_}, end_coordinate{end_coordinate_}, begin_coordinate{}, trace{}
    {}

    /*!\brief Constructor with identifier and score only.
     * \param id_               The alignment identifier.
     * \param score_            The alignment score.
     */
    align_result(id_t const id_,
                 score_t const score_)
        : id{id_}, score{score_}, end_coordinate{}, begin_coordinate{}, trace{}
    {}

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

    /*!\name Access functions
     * \brief Functions to access elements of the alignment result type.
     * \{
     */

    /*!\brief Returns the alignment identifier.
     * \return The id field.
     */
    constexpr id_t get_id() const noexcept
    {
        return id;
    }

    /*!\brief Returns the alignment score.
     * \return The score field.
     */
    constexpr score_t get_score() const noexcept
    {
        return score;
    }

    /*!\brief Returns the end coordinate of the alignment.
     * \return A pair of positions in the respective sequences, where the calculated alignment ends.
     * \throws std::bad_optional_access if the end coordinate calculation has not been requested.
     */
    constexpr coordinate_t get_end_coordinate() const
    {
        return end_coordinate.value();
    }

    /*!\brief Returns the begin coordinate of the alignment.
     * \return  A pair of positions in the respective sequences, where the calculated alignment starts.
     * \throws  std::bad_optional_access if the begin coordinate calculation has not been requested.
     * \details Guaranteed to be smaller than or equal to `get_end_coordinate()`.
     */
    constexpr coordinate_t get_begin_coordinate() const
    {
        return begin_coordinate.value();
    }

    /*!\brief Returns the traceback of the alignment.
     * \return At least two gapped sequences, which represent the alignment.
     * \throws std::bad_optional_access if the trace calculation has not been requested.
     */
    constexpr trace_t get_trace() const
    {
        return trace.value();
    }
    //!\}
};

} // namespace seqan3
