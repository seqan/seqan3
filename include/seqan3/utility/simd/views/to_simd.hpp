// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::detail::to_simd view.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <iterator>
#include <ranges>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/detail/template_inspection.hpp>
#include <seqan3/core/range/detail/adaptor_from_functor.hpp>
#include <seqan3/core/range/type_traits.hpp>
#include <seqan3/utility/simd/algorithm.hpp>
#include <seqan3/utility/simd/concept.hpp>
#include <seqan3/utility/simd/simd.hpp>
#include <seqan3/utility/simd/simd_traits.hpp>
#include <seqan3/utility/views/type_reduce.hpp>
#include <seqan3/utility/views/zip.hpp>

namespace seqan3::detail
{

/*!\brief Transforms a range of ranges into chunks of seqan3::simd vectors.
 * \implements std::ranges::input_range
 * \tparam urng_t The underlying range type; must model std::ranges::forward_range.
 * \tparam simd_t The simd type to convert to.
 * \ingroup utility_simd_views
 *
 * \details
 *
 * This view applies an Array-of-Structure to Structure-of-Array transformation on a set of sequences.
 * The typical application is to transform the memory layout of sequences such that they can be efficiently
 * used in vectorised algorithms. Accordingly, this view transforms the memory already into chunks of simd
 * vectors. The number of sequences must correspond to the length of the target simd vector, i.e. there are exactly as
 * many sequences as their are simd vector elements. The output range will be a range over chunks, where each
 * chunk represents a quadratic std::array over simd vectors (the size of the chunk is equal to the number of elements
 * in the target simd vector).
 *
 * Depending on the types of the input ranges a more efficient transformation using simd instructions is used.
 * The following requirements must be fulfilled by the inner range type of the underlying range:
 *  * they must model std::ranges::contiguous_range
 *  * the iterator and sentinel types must model std::sized_sentinel_for
 *  * the size of the rank type of the underlying alphabet must be 1.
 *
 * If one of these requirements is not fulfilled a standard fallback algorithm is used which might be slower to
 * transform the sequences, depending on the auto-vectorisation capabilities of the used compiler.
 */
template <std::ranges::view urng_t, simd::simd_concept simd_t>
class view_to_simd : public std::ranges::view_interface<view_to_simd<urng_t, simd_t>>
{
private:
    static_assert(std::ranges::forward_range<urng_t>, "The underlying range must model forward_range.");
    static_assert(std::ranges::input_range<std::ranges::range_value_t<urng_t>>,
                  "Expects the value type of the underlying range to be an input_range.");
    static_assert(std::default_initializable<std::ranges::iterator_t<std::ranges::range_value_t<urng_t>>>,
                  "Expects the inner range iterator to be default initializable.");
    static_assert(std::default_initializable<std::ranges::sentinel_t<std::ranges::range_value_t<urng_t>>>,
                  "Expects the inner range sentinel to be default initializable.");
    static_assert(semialphabet<std::ranges::range_value_t<std::ranges::range_value_t<urng_t>>>,
                  "Expects semi-alphabet as value type of the inner range.");

    /*!\name Auxiliary types
     * \{
     */
    using inner_range_type = std::ranges::range_reference_t<urng_t>;    //!< The inner range type.
    using chunk_type = std::array<simd_t, simd_traits<simd_t>::length>; //!< The underlying type to hold the chunks.
    using scalar_type = typename simd_traits<simd_t>::scalar_type;      //!< The scalar type.
    //!\brief The SIMD type with maximal number of lanes for the current arch.
    using max_simd_type = simd_type_t<uint8_t, simd_traits<simd_t>::max_length>;
    //!\}

    /*!\name Auxiliary variables
     * \{
     */
    //!\brief Check if fast load is enabled.
    static constexpr bool fast_load =
        std::ranges::contiguous_range<inner_range_type>
        && std::sized_sentinel_for<std::ranges::iterator_t<inner_range_type>, std::ranges::sentinel_t<inner_range_type>>
        && sizeof(alphabet_rank_t<std::ranges::range_value_t<inner_range_type>>) == 1;

    //!\brief The size of one chunk. Equals the number of elements in the simd vector.
    static constexpr uint8_t chunk_size = simd_traits<simd_t>::length;
    //!\brief The number of chunks that can be gathered with a single load.
    static constexpr uint8_t chunks_per_load = simd_traits<simd_t>::max_length / chunk_size;
    //!\brief The total number of chunks that can be cached.
    static constexpr uint8_t total_chunks = fast_load ? (chunks_per_load * chunks_per_load) : 1;
    //!\brief The alphabet size.
    static constexpr auto alphabet_size = seqan3::alphabet_size<std::ranges::range_value_t<inner_range_type>>;
    //!\}

    // Forward declare class' iterator type. See definition below.
    class iterator_type;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr view_to_simd()
        requires std::default_initializable<urng_t>
    = default;                                                          //!< Defaulted.
    constexpr view_to_simd(view_to_simd const &) = default;             //!< Defaulted.
    constexpr view_to_simd(view_to_simd &&) = default;                  //!< Defaulted.
    constexpr view_to_simd & operator=(view_to_simd const &) = default; //!< Defaulted.
    constexpr view_to_simd & operator=(view_to_simd &&) = default;      //!< Defaulted.
    ~view_to_simd() = default;                                          //!< Defaulted.

    /*!\brief Construction from the underlying range.
     * \param[in] urng The underlying range.
     * \param[in] padding_value The value used to fill up smaller sequences.
     */
    constexpr view_to_simd(urng_t urng, scalar_type const padding_value = alphabet_size) :
        urng{std::move(urng)},
        padding_simd_vector{simd::fill<simd_t>(padding_value)},
        padding_value{padding_value}
    {
        // Check if the size is less or equal the simd size.
        if (std::ranges::distance(urng) > chunk_size)
            throw std::invalid_argument{"The size of the underlying range must be less than or equal to the size of "
                                        "the given simd type!"};
    }

    //!\overload
    template <typename other_urng_t>
        requires (!std::same_as<std::remove_cvref_t<other_urng_t>, view_to_simd>)
              && (!std::same_as<other_urng_t, urng_t>) && std::ranges::viewable_range<other_urng_t>
    constexpr view_to_simd(other_urng_t && urng, scalar_type const padding_value = alphabet_size) :
        view_to_simd{views::type_reduce(std::forward<other_urng_t>(urng)), padding_value}
    {}
    //!\}

    /*!\name Iterators
     * \{
     */
    //!\brief The iterator to begin of this range.
    constexpr iterator_type begin() noexcept
    {
        return {*this};
    }

    //!\brief Const iteration is disabled.
    constexpr void begin() const noexcept = delete;

    //!\brief A sentinel representing the end of this range.
    constexpr std::default_sentinel_t end() noexcept
    {
        return std::default_sentinel;
    }

    //!\brief Const iteration is disabled.
    constexpr void end() const noexcept = delete;
    //!\}

    //!\brief Checks whether the range is empty.
    constexpr bool empty() const noexcept
        requires std::ranges::forward_range<inner_range_type>
    {
        return std::ranges::all_of(urng,
                                   [](auto & rng)
                                   {
                                       return std::ranges::empty(rng);
                                   });
    }

    /*!\brief Returns the size of this range.
     *
     * \details
     *
     * Only available if the inner range types model std::ranges::sized_range.
     */
    constexpr size_t size() const noexcept
        requires std::ranges::sized_range<inner_range_type>
    {
        auto it = std::ranges::max_element(urng,
                                           [](auto & lhs, auto & rhs)
                                           {
                                               return std::ranges::size(lhs) < std::ranges::size(rhs);
                                           });

        return (it != std::ranges::end(urng)) ? (std::ranges::size(*it) + chunk_size - 1) / chunk_size : 0;
    }

private:
    urng_t urng{};                                             //!< The underlying range.
    std::array<chunk_type, total_chunks> cached_simd_chunks{}; //!< The cached chunks of transformed simd vectors.
    simd_t padding_simd_vector{};                              //!< A cached simd vector with the padding symbol.
    scalar_type padding_value{}; //!< The padding value used to fill the corresponding simd vector element.
};

/*!\brief Iterator that transposes the underlying range of ranges and transforms that to SIMD types.
 *
 * \details
 *
 * This iterator models std::input_iterator.
 * When dereferencing it returns a span over a range with value type `simd_t`.
 */
template <std::ranges::view urng_t, simd::simd_concept simd_t>
class view_to_simd<urng_t, simd_t>::iterator_type
{
public:
    /*!\name Associated types
     * \{
     */
    using reference = std::span<std::ranges::range_value_t<chunk_type>>; //!< The reference type.
    using value_type = reference;                                        //!< The value type.
    using pointer = void;                                                //!< The pointer type.
    using difference_type = ptrdiff_t;                                   //!< The difference type.
    using iterator_category = std::input_iterator_tag;                   //!< The iterator category.
    using iterator_concept = iterator_category;                          //!< The iterator concept.
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr iterator_type() = default;                                  //!< Defaulted.
    constexpr iterator_type(iterator_type const &) = default;             //!< Defaulted.
    constexpr iterator_type(iterator_type &&) = default;                  //!< Defaulted.
    constexpr iterator_type & operator=(iterator_type const &) = default; //!< Defaulted.
    constexpr iterator_type & operator=(iterator_type &&) = default;      //!< Defaulted.
    ~iterator_type() = default;                                           //!< Defaulted.

    /*!\brief Construction from the associated range.
     * \param this_view A reference to the associated view.
     *
     * \details
     *
     * Initialises the iterators and sentinels of the underlying sequences to be transformed and calls
     * underflow to fetch the first chunk.
     */
    constexpr iterator_type(view_to_simd & this_view) : this_view{&this_view}, current_chunk_pos{0}
    {
        // Initialise the iterator of the sub ranges.
        size_t seq_id = 0;
        for (auto it = std::ranges::begin(this_view.urng); it != std::ranges::end(this_view.urng); ++it, ++seq_id)
        {
            cached_iter[seq_id] = std::ranges::begin(*it);
            cached_sentinel[seq_id] = std::ranges::end(*it);
        }

        // The batch is empty and by default the constructed iterator is pointing to the end.
        if (seq_id == 0)
            return;

        // The batch is not empty but might not be full either.
        // If a slot is supposed to be empty, it will be initialised with the iterator of the first sequence set to the
        // end emulating an empty sequence.
        auto sentinel_it = std::ranges::next(cached_iter[0], cached_sentinel[0]);
        for (; seq_id < chunk_size; ++seq_id)
        {
            cached_iter[seq_id] = sentinel_it;
            cached_sentinel[seq_id] = cached_sentinel[0];
        }

        // Check if this is the final chunk already.
        final_chunk = all_iterators_reached_sentinel();

        // Fetch the next available input characters from the sequences and transform them into simd vectors.
        underflow();
    }
    //!\}

    /*!\name Element access
     * \{
     */
    //!\brief Returns a reference to the current chunk of simd vectors.
    constexpr reference operator*() const noexcept
    {
        assert(this_view != nullptr);
        return std::span{this_view->cached_simd_chunks[current_chunk_pos].begin(),
                         (current_chunk_pos == final_chunk_pos) ? final_chunk_size : chunk_size};
    }
    //!\}

    /*!\name Arithmetic operators
     * \{
     */
    //!\brief Advances the iterator to the next chunk.
    constexpr iterator_type & operator++(/*pre-increment*/)
    {
        if constexpr (fast_load)
        { // Check if cached chunks have been already consumed and we need to fetch the next chunks.
            if (current_chunk_pos == final_chunk_pos)
            {
                underflow();
                current_chunk_pos = 0;
            }
            else
            {
                ++current_chunk_pos;
            }
        }
        else // In case fast load is not available only one chunk is filled at a time.
        {
            underflow();
        }

        return *this;
    }

    //!\brief Advances the iterator to the next chunk and returns the previous pointed-to value.
    constexpr value_type operator++(int /*post-increment*/)
    {
        value_type tmp = this->operator*();
        ++(*this);
        return tmp;
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */
    //!\brief Returns `true` if iterator reached the end, otherwise `false`.
    constexpr bool operator==(std::default_sentinel_t const &) const noexcept
    {
        return at_end;
    }

    //!\copydoc seqan3::detail::view_to_simd::iterator_type::operator==
    friend constexpr bool operator==(std::default_sentinel_t const &, iterator_type const & rhs) noexcept
    {
        return rhs.at_end;
    }

    //!\brief Returns `true` if iterator did not reach the end yet, otherwise `false`.
    constexpr bool operator!=(std::default_sentinel_t const &) const noexcept
    {
        return !at_end;
    }

    //!\copydoc seqan3::detail::view_to_simd::iterator_type::operator!=
    friend constexpr bool operator!=(std::default_sentinel_t const &, iterator_type const & rhs) noexcept
    {
        return !rhs.at_end;
    }
    //!\}

private:
    /*!\brief Unpacks one row of the transposed byte matrix using simd instructions.
     * \param[in] row The current matrix row to unpack.
     * \returns The unpacked simd vectors depending on the target vector length.
     *
     * \details
     *
     * If the target vector length is the same as the max vector length nothing will be unpacked and the
     * function is a noop. In the other cases the corresponding parts of the source vector are extracted and
     * upcasted to the target simd type.
     */
    auto unpack(max_simd_type const & row) const
    {
        if constexpr (chunk_size == simd_traits<max_simd_type>::length / 2) // upcast into 2 vectors.
        {
            return std::array{simd::upcast<simd_t>(extract_half<0>(row)),  // 1. half
                              simd::upcast<simd_t>(extract_half<1>(row))}; // 2. half
        }
        else if constexpr (chunk_size == simd_traits<max_simd_type>::length / 4) // upcast into 4 vectors.
        {
            return std::array{simd::upcast<simd_t>(extract_quarter<0>(row)),  // 1. quarter
                              simd::upcast<simd_t>(extract_quarter<1>(row)),  // 2. quarter
                              simd::upcast<simd_t>(extract_quarter<2>(row)),  // 3. quarter
                              simd::upcast<simd_t>(extract_quarter<3>(row))}; // 4. quarter
        }
        else if constexpr (chunk_size == simd_traits<max_simd_type>::length / 8) // upcast into 8 vectors.
        {
            return std::array{simd::upcast<simd_t>(extract_eighth<0>(row)),  // 1. eighth
                              simd::upcast<simd_t>(extract_eighth<1>(row)),  // 2. eighth
                              simd::upcast<simd_t>(extract_eighth<2>(row)),  // 3. eighth
                              simd::upcast<simd_t>(extract_eighth<3>(row)),  // 4. eighth
                              simd::upcast<simd_t>(extract_eighth<4>(row)),  // 5. eighth
                              simd::upcast<simd_t>(extract_eighth<5>(row)),  // 6. eighth
                              simd::upcast<simd_t>(extract_eighth<6>(row)),  // 7. eighth
                              simd::upcast<simd_t>(extract_eighth<7>(row))}; // 8. eighth
        }
        else
        {
            return std::array{simd::upcast<simd_t>(row)};
        }
    }

    /*!\brief Unpacks the matrix of simd types and caches the respective chunk entries.
     * \param matrix The transposed byte matrix from the efficient load procedure.
     *
     * \details
     *
     * In the efficient load procedure a quadratic byte matrix is first filled and then transposed using
     * efficient simd instructions. Depending on the target simd vector type this byte matrix must be unpacked
     * and the corresponding simd vectors must be assigned to their respective position within the cached
     * chunk array.
     */
    constexpr void split_into_sub_matrices(std::array<max_simd_type, simd_traits<max_simd_type>::length> matrix) const
    {
        auto apply_padding = [this](simd_t const vec)
        {
            return (vec == simd::fill<simd_t>(static_cast<uint8_t>(~0))) ? this_view->padding_simd_vector : vec;
        };

        // Iterate over the rows of the matrix
        for (uint8_t row = 0; row < static_cast<uint8_t>(matrix.size()); ++row)
        {
            // split a row into multiple chunks of size `chunk_size`
            auto chunked_row = unpack(matrix[row]);

            if constexpr (chunked_row.size() == 1)
            {
                this_view->cached_simd_chunks[0][row] = apply_padding(std::move(chunked_row[0]));
            }
            else // Parse the tuple elements and store them in the cached simd chunks.
            {
                static_assert(chunked_row.size() == chunks_per_load, "Expected chunks_per_load many simd vectors.");

                for (uint8_t chunk = 0; chunk < chunks_per_load; ++chunk) // store chunks in respective cached entries.
                {
                    size_t idx = chunk * chunks_per_load + row / chunk_size;
                    this_view->cached_simd_chunks[idx][row % chunk_size] = apply_padding(std::move(chunked_row[chunk]));
                }
            }
        }
    }

    /*!\brief Checks if all sequence iterators reached the end.
     * \returns `true` if all iterators reached the end, otherwise `false`.
     */
    constexpr bool all_iterators_reached_sentinel() const noexcept
    {
        using std::get;

        return std::ranges::all_of(views::zip(cached_iter, cached_sentinel),
                                   [](auto && iterator_sentinel_pair)
                                   {
                                       return get<0>(iterator_sentinel_pair) == get<1>(iterator_sentinel_pair);
                                   });
    }

    /*!\brief Convert a single column into a simd vector.
     * \tparam indices A non-type template parameter pack over the sequence indices.
     *
     * \returns A simd vector filled with the symbols of the sequences at the current position.
     *
     * \details
     *
     * Converts a single column over the sequences into a simd type. If the end of one sequence was already
     * reached it will return the padding value instead.
     */
    constexpr simd_t convert_single_column() noexcept
    {
        simd_t simd_column{};
        for (size_t idx = 0u; idx < chunk_size; ++idx)
        {
            if (cached_iter[idx] == cached_sentinel[idx])
            {
                simd_column[idx] = this_view->padding_value;
            }
            else
            {
                simd_column[idx] = static_cast<scalar_type>(seqan3::to_rank(*cached_iter[idx]));
                ++cached_iter[idx];
            }
        };
        return simd_column;
    }

    /*!\brief Updates the end of the final chunk and sets the index of the final chunk.
     *
     * \tparam array_t The array type containing the iterators over sequences.
     * \param[in] iterators_before_update The array containing the iterators before the load operations.
     *
     * \details
     *
     * Sets the index of the final chunk (`final_chunk_pos`) and updates the end position of the final chunk such that
     * the view ends at the last character of the longest sequences contained in the set of sequences to be transformed.
     */
    template <typename array_t>
    constexpr void update_final_chunk_position(array_t const & iterators_before_update) noexcept
    {
        size_t max_distance = 0;
        for (auto && [it, sent] : views::zip(iterators_before_update, cached_sentinel))
            max_distance = std::max<size_t>(std::ranges::distance(it, sent), max_distance);

        assert(max_distance > 0);
        assert(max_distance <= (total_chunks * chunk_size));

        --max_distance;
        final_chunk_pos = max_distance / chunk_size;
        // first we should be able to check the chunk position.
        final_chunk_size = (max_distance % chunk_size) + 1;
    }

    //!\brief Fetches the next available chunk(s).
    constexpr void underflow()
        requires fast_load
    {
        at_end = final_chunk;
        if (at_end) // reached end of stream.
            return;
        // For the efficient load we assume at most one byte sized alphabets.
        // Hence we can load `simd_traits<simd_t>::max_length` length many elements at once.
        // Depending on the packing of `simd_t` we can prefetch blocks and store them in the `cached_simd_chunks`.
        // E.g. assume `simd_t` with length 8 on SSE4 with max length 16.
        // To fill the 16x16 matrix we need four 8x8 matrices.
        // Thus, for the 8 sequences we need to load two times 16 consecutive bytes to fill the matrix, i.e. two loads
        // see figure below.
        //
        //       0    1  ...    7 |   8    9  ...   15
        // 0  [a00, a01, ..., a07]|[a08, a09, ..., a15]  // first load of seq a reads 16 characters
        // 1  [b00, b01, ..., b07]|[b08, b09, ..., b15]  // first load of seq b reads 16 characters
        //            ...         |        ...
        // 7  [g00, g01, ..., g07]|[g08, g09, ..., g15]  // first load of seq g reads 16 characters
        //    ----------------------------------------
        // 8  [a16, a17, ..., a23]|[a24, a25, ..., a31]  // second load of seq a reads next 16 characters
        // 9  [b16, b17, ..., b23]|[b24, b25, ..., b31]  // second load of seq b reads next 16 characters
        //            ...         |        ...
        // 15 [g16, g17, ..., g23]|[g24, g25, ..., g31]  // second load of seq g reads next 16 characters
        //
        // This quadratic byte matrix can be transposed efficiently with simd instructions.
        // If the target simd scalar type is bigger we can apply the same mechanism but have then 16 4x4 matrices
        // (32 bit) or 256 2x2 matrices (64 bit).

        constexpr int8_t max_size = simd_traits<simd_t>::max_length;
        std::array<max_simd_type, max_size> matrix{};
        decltype(cached_iter) iterators_before_update{cached_iter}; // Keep track of iterators before the update.
        // Iterate over each sequence.
        for (uint8_t sequence_pos = 0; sequence_pos < chunk_size; ++sequence_pos)
        { // Iterate over each block depending on the packing of the target simd vector.
            for (uint8_t chunk_pos = 0; chunk_pos < chunks_per_load; ++chunk_pos)
            {
                uint8_t pos = chunk_pos * chunk_size + sequence_pos; // matrix entry to fill
                if (cached_sentinel[sequence_pos] - cached_iter[sequence_pos] >= max_size)
                { // Not in final block, thus load directly from memory.
                    matrix[pos] = simd::load<max_simd_type>(std::addressof(*cached_iter[sequence_pos]));
                    std::advance(cached_iter[sequence_pos], max_size);
                }
                else // Loads the final block byte wise in order to not load from uninitialised memory.
                {
                    matrix[pos] = simd::fill<max_simd_type>(~0);
                    auto & sequence_it = cached_iter[sequence_pos];
                    for (int8_t idx = 0; sequence_it != cached_sentinel[sequence_pos]; ++sequence_it, ++idx)
                        matrix[pos][idx] = seqan3::to_rank(*sequence_it);
                }
            }
        }

        // Handle final chunk which might not end at an offset which is not a multiple of `chunk_size`.
        final_chunk = all_iterators_reached_sentinel();

        if (final_chunk)
            update_final_chunk_position(iterators_before_update);

        simd::transpose(matrix);
        split_into_sub_matrices(std::move(matrix));
    }

    //!\overload
    constexpr void underflow()
        requires (!fast_load)
    {
        at_end = final_chunk;
        if (at_end) // reached end of stream.
            return;

        decltype(cached_iter) iterators_before_update{cached_iter}; // Keep track of iterators before the update.
        for (size_t i = 0; i < chunk_size; ++i)
            this_view->cached_simd_chunks[0][i] = convert_single_column();

        final_chunk = all_iterators_reached_sentinel();

        if (final_chunk)
            update_final_chunk_position(iterators_before_update);
    }

    //!\brief Array containing the cached sequence iterators over the inner ranges.
    std::array<std::ranges::iterator_t<inner_range_type>, chunk_size> cached_iter{};
    //!\brief Array containing the cached sequence sentinels over the inner ranges.
    std::array<std::ranges::sentinel_t<inner_range_type>, chunk_size> cached_sentinel{};
    //!\brief Pointer to the associated range.
    view_to_simd * this_view{nullptr};
    //!\brief The size of the final chunk.
    uint8_t final_chunk_size{chunk_size};
    //!\brief The final chunk position.
    uint8_t final_chunk_pos{total_chunks - 1};
    //!\brief The current chunk position.
    uint8_t current_chunk_pos{0};
    //!\brief Flag indicating that final chunk was reached.
    bool final_chunk{true};
    //!\brief Flag indicating that iterator is at end.
    bool at_end{true};
};

// ============================================================================
//  to_simd_fn (adaptor definition)
// ============================================================================

/*!\brief views::to_simd's range adaptor closure object type.
 * \ingroup utility_simd_views
 * \tparam simd_t The target simd type.
 *
 * \details
 *
 * Returns a seqan3::detail::view_to_simd view for a given std::ranges::viewable_range.
 */
template <simd::simd_concept simd_t>
struct to_simd_fn
{
    //!\brief The type of a padding value.
    using padding_t = typename simd_traits<simd_t>::scalar_type;

    /*!\brief Returns a range adaptor closure object with the given parameter.
     * \param[in] padding_value The padding character to use for smaller sequences.
     */
    constexpr auto operator()(padding_t const padding_value) const noexcept
    {
        return detail::adaptor_from_functor{*this, padding_value};
    }

    //!\brief Returns a range adaptor object.
    constexpr auto operator()() const noexcept
    {
        return detail::adaptor_from_functor{*this};
    }

    /*!\brief Call the view's constructor with the underlying std::ranges::viewable_range as argument.
     * \param[in] urange The input range to process; must model std::ranges::forward_range and std::ranges::viewable_range.
     * \param[in] padding_value The padding character to use for smaller sequences.
     * \returns A range that transforms a collection of sequences into chunks of simd vectors.
     */
    template <std::ranges::range urng_t>
    constexpr auto operator()(urng_t && urange, padding_t const padding_value) const noexcept
    {
        static_assert(std::ranges::forward_range<urng_t>,
                      "The underlying range in views::to_simd must model std::ranges::forward_range.");
        static_assert(std::ranges::viewable_range<urng_t>,
                      "The underlying range in views::to_simd must model std::ranges::viewable_range.");
        static_assert(std::ranges::input_range<std::ranges::range_value_t<urng_t>>,
                      "The value type of the underlying range must model std::ranges::input_range.");
        static_assert(semialphabet<std::ranges::range_value_t<std::ranges::range_value_t<urng_t>>>,
                      "The value type of the inner ranges must model seqan3::semialphabet.");

        return view_to_simd<type_reduce_t<urng_t>, simd_t>{std::forward<urng_t>(urange), padding_value};
    }

    /*!\brief Call the view's constructor with the underlying std::ranges::viewable_range as argument.
     * \param[in] urange The input range to process; must model std::ranges::forward_range and std::ranges::viewable_range.
     * \returns A range that transforms a collection of sequences into chunks of simd vectors.
     */
    template <std::ranges::range urng_t>
    constexpr auto operator()(urng_t && urange) const noexcept
    {
        static_assert(std::ranges::forward_range<urng_t>,
                      "The underlying range in views::to_simd must model std::ranges::forward_range.");
        static_assert(std::ranges::viewable_range<urng_t>,
                      "The underlying range in views::to_simd must model std::ranges::viewable_range.");
        static_assert(std::ranges::input_range<std::ranges::range_value_t<urng_t>>,
                      "The value type of the underlying range must model std::ranges::input_range.");
        static_assert(semialphabet<std::ranges::range_value_t<std::ranges::range_value_t<urng_t>>>,
                      "The value type of the inner ranges must model seqan3::semialphabet.");

        return view_to_simd<type_reduce_t<urng_t>, simd_t>{std::forward<urng_t>(urange)};
    }

    //!\brief Overloaded bit-operator to allow chaining with other ranges.
    template <std::ranges::range urng_t>
    constexpr friend auto operator|(urng_t && urange, to_simd_fn const & me)
    {
        return me(std::forward<urng_t>(urange));
    }
};

} // namespace seqan3::detail

namespace seqan3::views
{

/*!\brief A view that transforms a range of ranges into chunks of seqan3::simd vectors.
 * \ingroup utility_simd_views
 * \tparam urng_t The type of the range being processed.
 * \tparam simd_t The target simd vector type.
 * \param[in] urange The range being processed.
 * \param[in] padding An optional padding value.
 * \returns A range of ranges with the original sequences transformed into simd vectors.
 *
 * \details
 *
 * \header_file{seqan3/utility/simd/views/to_simd.hpp}
 *
 * This view can be used to transform a collection of sequences into chunks of simd vectors. This transformation is
 * also known as Array-of-Structure to Structure-of-Array transformation. It is used to transform the memory layout of
 * the sequences to a more efficient form when used in vectorised algorithms. The number of sequences contained in the
 * range to be transformed cannot be larger than the number of elements stored in the target simd vector, i.e.
 * the size of `urange` <= simd_traits<simd_t>::length.
 * After applying the transformation one column of the outer range is transposed into a simd vector.
 * This means that the characters of all sequences at a given position `x` are stored in a simd vector retaining the
 * original order. The returned range itself is a range-of-ranges. When dereferencing the iterator a std::span over
 * a std::array with at most simd length many vectors is returned. If a sequence is empty or ends before the largest
 * sequence in the collection, it can be padded with an optional value.
 *
 * ### View properties
 *
 * | Concepts and traits              | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                                                |
 * |----------------------------------|:-------------------------------------:|:-----------------------------------------------------------------------------:|
 * | std::ranges::input_range         | *required*                            | *preserved*                                                                   |
 * | std::ranges::forward_range       | *required*                            | *lost*                                                                        |
 * | std::ranges::bidirectional_range |                                       | *lost*                                                                        |
 * | std::ranges::random_access_range |                                       | *lost*                                                                        |
 * | std::ranges::contiguous_range    |                                       | *lost*                                                                        |
 * |                                  |                                       |                                                                               |
 * | std::ranges::viewable_range      | *required*                            | *guaranteed*                                                                  |
 * | std::ranges::view                |                                       | *guaranteed*                                                                  |
 * | std::ranges::sized_range         |                                       | *preserved*  (iff `std::ranges::sized_range<std::ranges::range_value_t<urng_t>>` is `true`) |
 * | std::ranges::common_range        |                                       | *lost*                                                                        |
 * | std::ranges::output_range        |                                       | *lost*                                                                        |
 * | seqan3::const_iterable_range     |                                       | *lost*                                                                        |
 * |                                  |                                       |                                                                               |
 * | std::ranges::range_reference_t   |                                       | std::span<simd_t>                                                             |
 *
 *
 * * `urng_t` is the type of the range modified by this view (input).
 * * the expression `std::ranges::input_range<std::ranges::range_value_t<urng_t>` must evaluate to `true`
 * * the expression `std::default_initializable<std::ranges::range_value_t<urng_t>>` must evaluate to `true`
 * * the expression `semialphabet<std::ranges::range_value_t<std::ranges::range_value_t<urng_t>>>` must evaluate to `true`
 * * `rrng_type` is the type of the range returned by this view.
 * * for more details, see \ref views.
 *
 * ### Example
 *
 * \include test/snippet/utility/simd/views/to_simd.cpp
 *
 * The output is as follows:
 *
 * ```console
 * Chunk 0:
 * [0,0,2,2,0,0,2,8]
 * [1,2,0,2,1,1,3,8]
 * [2,3,1,0,2,2,0,8]
 * [3,2,3,3,3,0,1,8]
 * [0,0,0,1,0,3,2,8]
 * [1,2,2,2,1,1,2,8]
 * [2,1,1,0,2,2,0,8]
 * [3,3,0,1,0,0,3,8]
 *
 * Chunk 1:
 * [0,0,1,2,1,1,2,8]
 * [1,1,2,2,2,2,2,8]
 * [2,2,0,0,2,0,3,8]
 * [3,2,2,1,0,1,0,8]
 * [0,0,1,3,1,3,0,8]
 * [1,1,2,0,2,0,0,8]
 * [2,3,0,2,3,2,1,8]
 * [0,0,2,1,0,1,1,8]
 *
 * Chunk 2:
 * [3,2,0,8,1,2,2,8]
 * [1,1,3,8,2,0,1,8]
 * [2,3,1,8,0,1,0,8]
 * [8,0,2,8,2,8,1,8]
 * [8,1,8,8,1,8,0,8]
 * [8,2,8,8,2,8,3,8]
 * [8,0,8,8,0,8,8,8]
 * [8,1,8,8,2,8,8,8]
 *
 * Chunk 3:
 * [8,3,8,8,1,8,8,8]
 * [8,8,8,8,3,8,8,8]
 * [8,8,8,8,0,8,8,8]
 * [8,8,8,8,1,8,8,8]
 * [8,8,8,8,2,8,8,8]
 * [8,8,8,8,0,8,8,8]
 * [8,8,8,8,2,8,8,8]
 * [8,8,8,8,1,8,8,8]
 * ```
 *
 * \hideinitializer
 */

template <simd::simd_concept simd_t>
inline constexpr auto to_simd = detail::to_simd_fn<simd_t>{};

} // namespace seqan3::views
