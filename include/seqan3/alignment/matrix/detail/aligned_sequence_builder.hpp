// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
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
#include <seqan3/alignment/band/static_band.hpp>
#include <seqan3/alignment/matrix/detail/matrix_coordinate.hpp>
#include <seqan3/alignment/matrix/trace_directions.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/core/type_traits/concept.hpp>
#include <seqan3/core/type_traits/lazy.hpp>
#include <seqan3/range/decorator/gap_decorator.hpp>
#include <seqan3/range/views/convert.hpp>
#include <seqan3/range/views/slice.hpp>
#include <seqan3/range/views/type_reduce.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

/*!\brief A transformation trait that returns the correct aligned sequence type for a given sequence type.
 * \ingroup alignment_matrix
 * \implements seqan3::transformation_trait
 *
 * \tparam range_t The range type to be transformed.
 *
 * \details
 *
 * This transformation trait helps to define the correct type representing the aligned sequence.
 * The aligned sequence type is defined as follows:
 *
 * * a seqan3::gap_decorator over a slice of an instance of `range_t`, if it is constructible, otherwise
 * * a std::vector over the seqan3::gapped alphabet over the value type of `range_t`.
 *
 * Note some ranges cannot be decorated with the seqan3::gap_decorator, e.g. a non-random access range. In this case
 * the fallback is to always use std::vector and copy the content from the passed range to the vector.
 */
template <std::ranges::viewable_range range_t>
struct make_aligned_sequence_type
{
    // The following expressions are used to check if the sequence types can be used as template arguments for the
    // seqan3::gap_decorator. Ranges that do not model std::random_access_range for instance cannot be augmented with
    // the gap_decorator and need to be copied instead.

    //!\brief The resulting aligned sequence type.
    using type = lazy_conditional_t<is_class_template_declarable_with_v<gap_decorator,
                                                                        decltype(std::declval<range_t>()
                                                                               | views::slice(0, 1))>,
                                    lazy<gap_decorator, decltype(std::declval<range_t>() | views::slice(0, 1))>,
                                    lazy<std::vector, gapped<std::ranges::range_value_t<range_t>>>>;
};

/*!\brief Builds the alignment for a given pair of sequences and the respective trace.
 * \ingroup alignment_matrix
 * \tparam fst_sequence_t The first sequence of the pairwise alignment; must model std::ranges::viewable_range.
 * \tparam sec_sequence_t The first sequence of the pairwise alignment; must model std::ranges::viewable_range.
 *
 * \details
 *
 * This class builds the alignment from a given trace path over the specified sequences. Use the interface
 * seqan3::detail::aligned_sequence_builder::operator() to get the actual alignment.
 * The returned seqan3::aligned_sequence type is determined by the input types `fst_sequence_t` and `sec_sequence_t`.
 *
 * See the seqan3::detail::make_aligned_sequence_type transformation trait for more information about the selected
 * type.
 *
 * Depending on the used alignment algorithm the computed alignment might only cover a subrange over the original
 * sequences. Accordingly, the returned alignment covers only the part of the sequences that are part of the given
 * trace path. One can use the seqan3::detail::aligned_sequence_builder::result_type to access the build alignment
 * and also to access the actual slice positions over which the alignment was built for the first sequence and
 * respectively the second sequence.
 */
template <std::ranges::viewable_range fst_sequence_t, std::ranges::viewable_range sec_sequence_t>
class aligned_sequence_builder
{
private:
    //!\brief The aligned sequence type for the first sequence.
    using fst_aligned_t = typename make_aligned_sequence_type<fst_sequence_t>::type;
    //!\brief The aligned sequence type for the second sequence.
    using sec_aligned_t = typename make_aligned_sequence_type<sec_sequence_t>::type;

    static_assert(seqan3::aligned_sequence<fst_aligned_t>,
                  "fst_aligned_t is required to model seqan3::aligned_sequence!");
    static_assert(seqan3::aligned_sequence<sec_aligned_t>,
                  "sec_aligned_t is required to model seqan3::aligned_sequence!");

public:

    //!\brief The result type when building the aligned sequences.
    struct [[nodiscard]] result_type
    {
        //!\brief The slice positions of the first sequence.
        std::pair<size_t, size_t> first_sequence_slice_positions{};
        //!\brief The slice positions of the second sequence.
        std::pair<size_t, size_t> second_sequence_slice_positions{};
        //!\brief The alignment over the slices of the first and second sequence, which corresponds to the given
        //!\      trace path.
        std::pair<fst_aligned_t, sec_aligned_t> alignment{};
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

    /*!\brief Construction from the underlying sequences.
     * \param[in] fst_rng The first range to build the aligned sequence for.
     * \param[in] sec_rng The second range to build the aligned sequence for.
     */
    constexpr aligned_sequence_builder(fst_sequence_t fst_rng, sec_sequence_t sec_rng) :
        fst_rng{views::type_reduce(std::forward<fst_sequence_t>(fst_rng))},
        sec_rng{views::type_reduce(std::forward<sec_sequence_t>(sec_rng))}
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

        result_type res{};
        auto trace_it = std::ranges::begin(trace_path);
        std::tie(res.first_sequence_slice_positions.second, res.second_sequence_slice_positions.second) =
            std::pair<size_t, size_t>{trace_it.coordinate()};

        std::vector<std::pair<trace_directions, size_t>> trace_segments;

        while (trace_it != std::ranges::end(trace_path))
        {
            trace_directions last_dir = *trace_it;
            size_t span = 0;
            for (; trace_it != std::ranges::end(trace_path) && *trace_it == last_dir; ++trace_it, ++span)
            {}

            trace_segments.emplace_back(last_dir, span);
        }

        std::tie(res.first_sequence_slice_positions.first, res.second_sequence_slice_positions.first) =
            std::pair<size_t, size_t>{trace_it.coordinate()};

        assign_unaligned(res.alignment.first, fst_rng  | views::slice(res.first_sequence_slice_positions.first,
                                                                      res.first_sequence_slice_positions.second));
        assign_unaligned(res.alignment.second, sec_rng | views::slice(res.second_sequence_slice_positions.first,
                                                                      res.second_sequence_slice_positions.second));

        // Now we need to insert the values.
        fill_aligned_sequence(trace_segments | std::views::reverse, res.alignment.first, res.alignment.second);

        return res;
    }

private:

    /*!\brief Fills the sequences with gaps according to the given trace segments.
     * \tparam reverse_traces_t The type storing the reverse trace.
     * \param[in] rev_traces The trace segments in order from source to sink in the trace matrix.
     * \param[in,out] fst_aligned The first aligned sequence to insert gaps into.
     * \param[in,out] sec_aligned The second aligned sequence to insert gaps into.
     */
    template <typename reverse_traces_t, typename fst_aligned_t, typename sec_aligned_t>
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
            if (dir == trace_directions::up)
                fst_it = insert_gap(fst_aligned, fst_it, span);

            if (dir == trace_directions::left)
                sec_it = insert_gap(sec_aligned, sec_it, span);

            fst_it += span;
            sec_it += span;
        }
    }

    type_reduce_view<fst_sequence_t> fst_rng; //!< A view over the first range.
    type_reduce_view<sec_sequence_t> sec_rng; //!< A view over the second range.
};

/*!\name Type deduction guides
 * \relates seqan3::detail::aligned_sequence_builder
 * \{
 */
//!\brief Deduces the type from the passed constructor arguments.
template <std::ranges::viewable_range fst_sequence_t, std::ranges::viewable_range sec_sequence_t>
aligned_sequence_builder(fst_sequence_t &&, sec_sequence_t &&) ->
    aligned_sequence_builder<fst_sequence_t, sec_sequence_t>;
//!\}
} // namespace seqan3::detail
