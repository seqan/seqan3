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

#include <optional>
#include <type_traits>
#include <vector>

#include <seqan3/alignment/aligned_sequence/aligned_sequence_concept.hpp>
#include <seqan3/alignment/band/static_band.hpp>
#include <seqan3/alignment/matrix/detail/matrix_coordinate.hpp>
#include <seqan3/alignment/matrix/trace_directions.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/core/type_traits/lazy.hpp>
#include <seqan3/range/decorator/gap_decorator.hpp>
#include <seqan3/range/views/slice.hpp>
#include <seqan3/range/views/view_all.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

namespace seqan3::detail
{
/*!\brief Builds the alignment for a given pair of sequences and the respective trace.
 * \ingroup alignment_matrix
 * \tparam fst_rng_t The first range of the pairwise alignment; must model std::ranges::viewable_range.
 * \tparam sec_rng_t The first range of the pairwise alignment; must model std::ranges::viewable_range.
 *
 * \details
 *
 * This class builds the alignment from a given trace path over the specified sequences. Use the interface
 * seqan3::detail::aligned_sequence_builder::operator() to get the alignment. The returned seqan3::aligned_sequence type
 * is determined by the input types `fst_rng_t` and `sec_rng_t`. If the respective template parameter can be wrapped
 * in a seqan3::gap_decorator after applying views::slice this type will be used, otherwise a std::vector speciallised
 * with seqan3::gapped alphabet over the underlying alphabet of the range will be used.
 */
template <std::ranges::viewable_range fst_rng_t, std::ranges::viewable_range sec_rng_t>
class aligned_sequence_builder
{
private:
    //!\brief The aligned sequence type for the first range.
    using fst_aligned_t =
        lazy_conditional_t<
            is_instantiable_with_v<gap_decorator, decltype(std::declval<fst_rng_t>() | views::slice(0, 1))>,
            lazy<gap_decorator, decltype(std::declval<fst_rng_t>() | views::slice(0, 1))>,
            lazy<std::vector, gapped<std::ranges::range_value_t<fst_rng_t>>>>;
    //!\brief The aligned sequence type for the second range.
    using sec_aligned_t =
        lazy_conditional_t<
            is_instantiable_with_v<gap_decorator, decltype(std::declval<sec_rng_t>() | views::slice(0, 1))>,
            lazy<gap_decorator, decltype(std::declval<sec_rng_t>() | views::slice(0, 1))>,
            lazy<std::vector, gapped<std::ranges::range_value_t<sec_rng_t>>>>;

    static_assert(seqan3::aligned_sequence<fst_aligned_t>,
                  "fst_aligned_t is required to model seqan3::aligned_sequence!");
    static_assert(seqan3::aligned_sequence<sec_aligned_t>,
                  "sec_aligned_t is required to model seqan3::aligned_sequence!");

public:

    //!\brief The result type when building the aligned sequences.
    struct result_type
    {
        matrix_coordinate begin{}; //!< The coordinate where the trimmed first and second sequence begin.
        matrix_coordinate end{}; //!< The coordinate where the trimmed first and second sequence end.
        std::pair<fst_aligned_t, sec_aligned_t> alignment{}; //!< The alignment over the trimmed sequences.
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
        fst_rng{views::all(std::forward<fst_rng_t>(fst_rng))},
        sec_rng{views::all(std::forward<sec_rng_t>(sec_rng))}
    {}
    //!\}

    /*!\brief Sets the column offset corresponding to the upper diagonal of the band within the matrix.
     * \param[in] offset The offset to set.
     *
     * \details
     *
     * The offset corresponds to the cell where the upper diagonal intersects with the first row of the matrix after
     * trimming the sequences. This option is only necessary if the matrix was computed with a banded alignment.
     * In this case, the row index of the alignment coordinate must be synchronized with the actual global row index
     * using the given offset.
     */
    constexpr void set_band_column_offset(size_t const offset)
    {
        band_col_offset = offset;
    }

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
        res.end = trace_it.coordinate();

        std::vector<std::pair<trace_directions, size_t>> traces;

        while (trace_it != std::ranges::end(trace_path))
        {
            trace_directions last_dir = *trace_it;
            size_t span = 0;
            for (; trace_it != std::ranges::end(trace_path) && *trace_it == last_dir; ++trace_it, ++span)
            {}

            traces.emplace_back(last_dir, span);
        }

        res.begin = trace_it.coordinate();

        if (band_col_offset.has_value())
        {
            res.begin.row += res.begin.col - band_col_offset.value();
            res.end.row += res.end.col - band_col_offset.value();
        }

        assign_unaligned(res.alignment.first, fst_rng | views::slice(res.begin.col, res.end.col));
        assign_unaligned(res.alignment.second, sec_rng | views::slice(res.begin.row, res.end.row));

        // Now we need to insert the values.
        fill_aligned_sequence(traces | std::views::reverse, res.alignment.first, res.alignment.second);

        return res;
    }

private:

    /*!\brief Fills the sequences with gaps according to the given trace segments.
     * \tparam reverse_traces_t The type storing the reverse trace.
     * \param[in] rev_traces The trace segments in order from source to sink in the trace matrix.
     * \param[in,out] fst_aligned The first aligned sequence to insert gaps into.
     * \param[in,out] sec_aligned The second aligned sequence to insert gaps into.
     */
    template <typename reverse_traces_t, typename fst_t, typename sec_t>
    void fill_aligned_sequence(reverse_traces_t && rev_traces, fst_t && fst_aligned, sec_t && sec_aligned) const
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

    all_view<fst_rng_t> fst_rng; //!< A view over the first range.
    all_view<sec_rng_t> sec_rng; //!< A view over the second range.
    std::optional<size_t> band_col_offset; //!< The offset of the column where an optional band starts.
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
