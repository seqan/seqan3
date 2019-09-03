// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::aligned_sequence_builder.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <type_traits>
#include <vector>

#include <seqan3/alignment/aligned_sequence/aligned_sequence_concept.hpp>
#include <seqan3/alignment/matrix/detail/matrix_coordinate.hpp>
#include <seqan3/alignment/matrix/trace_directions.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/range/decorator/gap_decorator.hpp>
#include <seqan3/range/view/slice.hpp>
#include <seqan3/range/view/view_all.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

/*!\brief Transforms the given range types into their aligned sequence format.
 * \ingroup alignment_matrix
 * \tparam fst_rng_t The type of the first range to be transformed to a aligned sequence type;
 *                   must model std::ranges::viewable_range.
 * \tparam sec_rng_t The type of the second range to be transformed to a aligned sequence type;
 *                   must model std::ranges::viewable_range.
 *
 * \details
 *
 * First this trait refines the original range types with seqan3::view::slice. Then it transforms the given range type
 * into a seqan3::gap_decorator if the expression is not ill-formed, otherwise the transformed aligned sequence type
 * is a std::vector over the respective alphabet decorated with seqan3::gapped.
 */
template <std::ranges::viewable_range fst_rng_t, std::ranges::viewable_range sec_rng_t>
struct aligned_sequence_builder_default_traits
{
private:
    //!\brief Helper type alias for sliced range.
    using fst_slice_t = decltype(std::declval<fst_rng_t>() | view::slice(0, 1));
    //!\copydoc aligned_sequence_builder_default_traits::fst_slice_t
    using sec_slice_t = decltype(std::declval<sec_rng_t>() | view::slice(0, 1));

    //!\brief Helper function to determine the transformed type.
    template <typename rng_t>
    //!\cond
        requires std::constructible_from<gap_decorator<rng_t>>
    //!\endcond
    static auto aligned_type() { return gap_decorator<rng_t>{}; }

    //!\overload
    template <typename rng_t>
    static auto aligned_type() { return std::vector<gapped<value_type_t<rng_t>>>{}; }

public:
    //!\brief The aligned sequence type for the first sequence.
    using fst_aligned_t = decltype(aligned_type<fst_slice_t>());
    //!\brief The aligned sequence type for the second sequence.
    using sec_aligned_t = decltype(aligned_type<sec_slice_t>());
};

/*!\brief Builds the alignment for a given pair of sequences and the respective trace.
 * \ingroup alignment_matrix
 * \tparam fst_rng_t The first range of the pairwise alignment; must model std::ranges::viewable_range.
 * \tparam sec_rng_t The first range of the pairwise alignment; must model std::ranges::viewable_range.
 * \tparam traits_t A traits class that determines the seqan3::aligned_sequence type; defaults to
 *                  seqan3::detail::aligned_sequence_builder_default_traits.
 *
 * \details
 *
 * This class builds the alignment from a given trace path over the specified sequences. Use the interface
 * seqan3::detail::aligned_sequence_builder::operator() to get the alignment. The returned seqan3::aligned_sequence type
 * is determined by the `traits_t` class which needs two member types that represent the seqan::aligned_sequence type.
 * The specified types must be constructible from the original type refined with a seqan3::view::slice. So the following
 * expressions must evaluate to `true` in order use the sequences with the given aligned sequence types:
 * ```cpp
 * std::constructible_from<typename traits_t::fst_aligned_t, decltype(std::declval<fst_rng_t>() | view::slice(0, 1))>
 * std::constructible_from<typename traits_t::sec_aligned_t, decltype(std::declval<sec_rng_t>() | view::slice(0, 1))>
 * ```
 */
template <std::ranges::viewable_range fst_rng_t,
          std::ranges::viewable_range sec_rng_t,
          typename traits_t = aligned_sequence_builder_default_traits<fst_rng_t, sec_rng_t>>
class aligned_sequence_builder
{
private:

    static_assert(std::constructible_from<typename traits_t::fst_aligned_t,
                                     decltype(std::declval<fst_rng_t>() | view::slice(0, 1))>,
                  "The given sequence type cannot be transformed to a aligned_sequence. It needs to be a view::slice "
                  "over the original type.");

    static_assert(std::constructible_from<typename traits_t::sec_aligned_t,
                                     decltype(std::declval<sec_rng_t>() | view::slice(0, 1))>,
                  "The given sequence type cannot be transformed to a aligned_sequence. It needs to be a view::slice "
                  "over the original type.");

    using fst_aligned_t = typename traits_t::fst_aligned_t; //!< The aligned sequence type for the first range.
    using sec_aligned_t = typename traits_t::sec_aligned_t; //!< The aligned sequence type for the second range.

public:

    //!\brief The result type when building the aligned sequences.
    struct result_type
    {
        matrix_coordinate begin; //!< The coordinate where the trimmed first and second sequence begin.
        matrix_coordinate end; //!< The coordinate where the trimmed first and second sequence end.
        std::pair<fst_aligned_t, sec_aligned_t> alignment; //!< The alignment over the trimmed sequences.
    };

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr aligned_sequence_builder() = default; //!< Defaulted.
    constexpr aligned_sequence_builder(aligned_sequence_builder const &) = default; //!< Defaulted.
    constexpr aligned_sequence_builder(aligned_sequence_builder &&) = default; //!< Defaulted.
    constexpr aligned_sequence_builder & operator=(aligned_sequence_builder const &) = default; //!< Defaulted.
    constexpr aligned_sequence_builder & operator=(aligned_sequence_builder &&) = default; //!< Defaulted.
    ~aligned_sequence_builder() = default; //!< Defaulted.

    /*\brief Construction from the underlying sequences.
     * \param[in] fst_rng The first range to build the aligned sequence for.
     * \param[in] sec_rng The second range to build the aligned sequence for.
     */
    constexpr aligned_sequence_builder(fst_rng_t fst_rng, sec_rng_t sec_rng) :
        fst_rng{view::all(std::forward<fst_rng_t>(fst_rng))},
        sec_rng{view::all(std::forward<sec_rng_t>(sec_rng))}
    {}
    //!\}

    /*!\brief Builds the aligned sequences from the given trace path.
     * \tparam trace_path_t The type of the trace path; must model std::ranges::input_range and
     *                      std::same_as<value_type_t<trace_path_t>, seqan::detail::trace_directions> must evaluate to
     *                      `true`.
     * \param[in] trace_path The trace path.
     * \returns seqan3::detail::aligned_sequence_builder::result_type with the built alignment.
     *
     * \details
     *
     * From the given trace path this function builds the aligned sequences for the first and the second
     * target sequence. The return type seqan3::detail::aligned_sequence_builder::result_type is an aggregate
     * type containing the begin and end coordinates with the respective sequence positions and the respective alignment
     * over the sliced sequences.
     */
    template <std::ranges::input_range trace_path_t>
    result_type operator()(trace_path_t && trace_path)
    {
        static_assert(std::same_as<value_type_t<trace_path_t>, trace_directions>,
                      "The value type of the trace path must be seqan3::detail::trace_directions");

        auto trace_it = std::ranges::begin(trace_path);
        matrix_coordinate path_end = trace_it.coordinate();

        std::vector<std::pair<trace_directions, size_t>> traces;

        while (trace_it != std::ranges::end(trace_path))
        {
            trace_directions last_dir = *trace_it;
            size_t span = 0;
            for (; trace_it != std::ranges::end(trace_path) && *trace_it == last_dir; ++trace_it, ++span)
            {}

            traces.emplace_back(last_dir, span);
        }

        matrix_coordinate path_begin = trace_it.coordinate();

        fst_aligned_t fst_aligned{fst_rng | view::slice(path_begin.col, path_end.col)};
        sec_aligned_t sec_aligned{sec_rng | view::slice(path_begin.row, path_end.row)};

        // Now we need to insert the values.
        fill_aligned_sequence(traces | std::view::reverse, fst_aligned, sec_aligned);

        return {path_begin, path_end, std::pair{fst_aligned, sec_aligned}};
    }

private:

    /*!\brief Fills the sequences with gaps according to the given trace segments.
     * \tparam reverse_traces_t The type storing the reverse trace.
     * \param[in] rev_traces The trace segments in order from source to sink in the trace matrix.
     * \param[in,out] fst_aligned The first aligned sequence to insert gaps into.
     * \param[in,out] sec_aligned The second aligned sequence to insert gaps into.
     */
    template <typename reverse_traces_t>
    void fill_aligned_sequence(reverse_traces_t && rev_traces,
                               fst_aligned_t & fst_aligned,
                               sec_aligned_t & sec_aligned) const
    {
        if (std::ranges::empty(rev_traces))
            return;

        auto fst_it = std::ranges::begin(fst_aligned);
        auto sec_it = std::ranges::begin(sec_aligned);

        for (auto const & [dir, span] : rev_traces)
        {
            switch (dir)
            {
                case trace_directions::diagonal:
                {
                    std::ranges::advance(fst_it, span);
                    std::ranges::advance(sec_it, span);
                    break;
                }
                case trace_directions::up:
                {
                    fst_it = insert_gap(fst_aligned, fst_it, span);
                    std::ranges::advance(sec_it, span);
                    break;
                }
                case trace_directions::left:
                {
                    sec_it = insert_gap(sec_aligned, sec_it, span);
                    std::ranges::advance(fst_it, span);
                    break;
                }
                case trace_directions::left_open:
                case trace_directions::up_open:
                case trace_directions::none:
                default:
                {
                    assert(false); // This should never happen.
                }
            }
        }
    }

    all_view<fst_rng_t> fst_rng; //!< A view over the first range.
    all_view<sec_rng_t> sec_rng; //!< A view over the second range.
};

/*!\name Type deduction guides
 * \relates seqan3::detail::aligned_sequence_builder
 * \{
 */
//!\brief Deduces the type from the passed constructor arguments.
template <std::ranges::viewable_range fst_rng_t, std::ranges::viewable_range sec_rng_t>
aligned_sequence_builder(fst_rng_t, sec_rng_t) -> aligned_sequence_builder<fst_rng_t, sec_rng_t>;
//!\}
} // namespace seqan3::detail
