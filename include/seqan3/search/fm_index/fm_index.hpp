// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \brief Provides the unidirectional seqan3::fm_index.
 */

#pragma once

#include <sdsl/suffix_trees.hpp>

#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/range/views/join.hpp>
#include <seqan3/range/views/to_rank.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/search/fm_index/concept.hpp>
#include <seqan3/search/fm_index/detail/fm_index_cursor.hpp>
#include <seqan3/search/fm_index/fm_index_cursor.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>

namespace seqan3
{

//!\cond
// forward declarations
template <typename index_t>
class fm_index_cursor;

template <typename index_t>
class bi_fm_index_cursor;
//!\endcond

/*!\addtogroup submodule_fm_index
 * \{
 */

/*!\brief The FM Index Configuration using a Wavelet Tree.
 *
 * \details
 *
 * ### Running time / Space consumption
 *
 * \f$SAMPLING\_RATE = 16\f$ \n
 * \f$\Sigma\f$: alphabet_size<alphabet_type> where alphabet_type is the seqan3 alphabet type (e.g. seqan3::dna4 has an
 *               alphabet size of 4).
 *
 * For an index over a text collection a delimiter is added in between the texts. This causes sigma to increase by 1.
 * \attention For any alphabet, the symbol with rank 255 is not allowed to occur in the text. Additionally,
 *            rank 254 cannot occur when indexing text collections.
 *
 * \if DEV
 * Rank 255 is generally not allowed. When constructing the index, we increase the rank of every letter by 1, since the
 * SDSL uses 0 as a sentinel. Incrementing 255 by 1 will cause the rank (uint8_t) to overflow to 0.
 * For text collections, 255 is reserved as delimiter between the individual texts. Therefore the letter with rank 254
 * cannot occur in the text.
 *
 * This index will only work for byte alphabets, i.e. alphabets with a size <= 256.
 * When switching to bigger alphabets, the static_asserts, the delimiter choice and the transform view for the
 * construction need to be adjusted. Additionally, all occurrences of uint8_t should be double checked to make sure
 * they also apply to bigger alphabets.
 * \endif
 *
 * \f$T_{BACKWARD\_SEARCH}: O(\log \Sigma)\f$
 *
 * \if DEV \todo Asymptotic space consumption: \endif
 *
 */
using sdsl_wt_index_type =
    sdsl::csa_wt<sdsl::wt_blcd<sdsl::bit_vector,
                               sdsl::rank_support_v<>,
                               sdsl::select_support_scan<>,
                               sdsl::select_support_scan<0>>,
                 16,
                 10'000'000,
                 sdsl::sa_order_sa_sampling<>,
                 sdsl::isa_sampling<>,
                 sdsl::plain_byte_alphabet>;

/*!\brief The default FM Index Configuration.
 * \attention The default might be changed in a future release. If you rely on a stable API and on-disk-format,
 *            please hard-code your sdsl_index_type to a concrete type.
 */
using default_sdsl_index_type = sdsl_wt_index_type;

/*!\brief The SeqAn FM Index.
 * \implements seqan3::fm_index_specialisation
 * \tparam alphabet_t        The alphabet type; must model seqan3::semialphabet.
 * \tparam text_layout_mode_ Indicates whether this index works on a text collection or a single text.
 *                           See seqan3::text_layout.
 * \tparam sdsl_index_type_  The type of the underlying SDSL index, must model seqan3::sdsl_index.
 * \implements seqan3::cerealisable
 * \details
 *
 * The seqan3::fm_index is a fast and space-efficient string index to search strings and collections of strings.
 *
 * ### General information
 *
 * Here is a short example on how to build an index and search a pattern using an cursor. Please note that there is a
 * very powerful search module with a high-level interface seqan3::search that encapsulates the use of cursors.
 *
 * \include test/snippet/search/fm_index.cpp
 *
 * \attention When building an index for a **single text** over any alphabet, the symbol with rank 255 is reserved
 *            and may not occur in the text.
 *
 * Here is an example using a collection of strings (e.g. a genome with multiple chromosomes or a protein database):
 *
 * \include test/snippet/search/fm_index_collection.cpp
 *
 * \attention When building an index for a **text collection** over any alphabet, the symbols with rank 254 and 255
 *            are reserved and may not be used in the text.
 *
 * \if DEV
 * ### Choosing an index implementation
 *
 * The underlying implementation of the FM Index (rank data structure, sampling rates, etc.) can be specified by
 * passing a new SDSL index type as second template parameter:
 *
 * \todo Link to SDSL documentation or write our own once SDSL3 documentation is available somewhere....
 *
 * \endif
 */
template <semialphabet alphabet_t,
          text_layout text_layout_mode_,
          detail::sdsl_index sdsl_index_type_ = default_sdsl_index_type>
class fm_index
{
private:
    /*!\name Member types
     * \{
     */
    //!\brief The type of the underlying SDSL index.
    using sdsl_index_type = sdsl_index_type_;
    /*!\brief The type of the reduced alphabet type. (The reduced alphabet might be smaller than the original alphabet
     *        in case not all possible characters occur in the indexed text.)
     */
    using sdsl_char_type = typename sdsl_index_type::alphabet_type::char_type;
    //!\brief The type of the alphabet size of the underlying SDSL index.
    using sdsl_sigma_type = typename sdsl_index_type::alphabet_type::sigma_type;
    //!\}

    //!\brief Underlying index from the SDSL.
    sdsl_index_type index;

    //!\brief Bitvector storing begin positions for collections.
    sdsl::sd_vector<> text_begin;
    //!\brief Select support for text_begin.
    sdsl::select_support_sd<1> text_begin_ss;
    //!\brief Rank support for text_begin.
    sdsl::rank_support_sd<1> text_begin_rs;

    /*!\brief Constructs the index given a range.
              The range cannot be an rvalue (i.e. a temporary object) and has to be non-empty.
     * \tparam text_t The type of range to construct from; must model std::ranges::bidirectional_range.
     * \param[in] text The text to construct from.
     *
     * \details
     * \if DEV
     * \todo This has to be better implemented with regard to the memory peak due to not matching interfaces
     *       with the SDSL.
     * \endif
     *
     * ### Complexity
     *
     * \if DEV \todo \endif At least linear.
     *
     * ### Exceptions
     *
     * No guarantee. \if DEV \todo Ensure strong exception guarantee. \endif
     */
    template <std::ranges::range text_t>
    //!\cond
        requires text_layout_mode_ == text_layout::single
    //!\endcond
    void construct(text_t && text)
    {
        static_assert(std::ranges::bidirectional_range<text_t>, "The text must model bidirectional_range.");
        static_assert(alphabet_size<innermost_value_type_t<text_t>> <= 256, "The alphabet is too big.");
        static_assert(std::convertible_to<innermost_value_type_t<text_t>, alphabet_t>,
                     "The alphabet of the text collection must be convertible to the alphabet of the index.");
        static_assert(dimension_v<text_t> == 1, "The input cannot be a text collection.");

        // text must not be empty
        if (std::ranges::empty(text))
            throw std::invalid_argument("The text that is indexed cannot be empty.");

        constexpr auto sigma = alphabet_size<alphabet_t>;

        // TODO:
        // * check what happens in sdsl when constructed twice!
        // * choose between in-memory/external and construction algorithms
        // * sdsl construction currently only works for int_vector, std::string and char *, not ranges in general
        // uint8_t largest_char = 0;
        sdsl::int_vector<8> tmp_text(std::ranges::distance(text));

        std::ranges::move(text
                          | views::to_rank
                          | std::views::transform([] (uint8_t const r)
                          {
                              if constexpr (sigma == 256)
                              {
                                  if (r == 255)
                                      throw std::out_of_range("The input text cannot be indexed, because for full"
                                                              "character alphabets the last one/two values are reserved"
                                                              "(single sequence/collection).");
                              }
                              return r + 1;
                          })
                          | std::views::reverse,
                          std::ranges::begin(tmp_text)); // reverse and increase rank by one

        sdsl::construct_im(index, tmp_text, 0);

        // TODO: would be nice but doesn't work since it's private and the public member references are const
        // index.m_C.resize(largest_char);
        // index.m_C.shrink_to_fit();
        // index.m_sigma = largest_char;
    }

    //!\overload
    template <std::ranges::range text_t>
    //!\cond
        requires text_layout_mode_ == text_layout::collection
    //!\endcond
    void construct(text_t && text)
    {
        static_assert(std::ranges::bidirectional_range<text_t>, "The text collection must model bidirectional_range.");
        static_assert(std::ranges::bidirectional_range<std::ranges::range_reference_t<text_t>>,
                      "The elements of the text collection must model bidirectional_range.");
        static_assert(alphabet_size<innermost_value_type_t<text_t>> <= 256, "The alphabet is too big.");
        static_assert(std::convertible_to<innermost_value_type_t<text_t>, alphabet_t>,
                     "The alphabet of the text collection must be convertible to the alphabet of the index.");
        static_assert(dimension_v<text_t> == 2, "The input must be a text collection.");

        // text collection must not be empty
        if (std::ranges::begin(text) == std::ranges::end(text))
            throw std::invalid_argument("The text that is indexed cannot be empty.");

        size_t text_size{0}; // text size including delimiters

        // there must be at least one non-empty text in the collection
        bool all_empty = true;

        for (auto && t : text)
        {
            if (std::ranges::begin(t) != std::ranges::end(t))
            {
                all_empty = false;
            }
            text_size += 1 + std::ranges::distance(t); // text size and delimiter (sum will be 1 for empty texts)
        }

        if (all_empty)
            throw std::invalid_argument("A text collection that only contains empty texts cannot be indexed.");

        constexpr auto sigma = alphabet_size<alphabet_t>;

        // Instead of creating a bitvector of size `text_size`, setting the bits to 1 and then compressing it, we can
        // use the `sd_vector_builder(text_size, number_of_ones)` because we know the parameters and the 1s we want to
        // set are in a strictly increasing order. This inplace construction of the compressed vector saves memory.
        sdsl::sd_vector_builder builder(text_size, std::ranges::distance(text));
        size_t prefix_sum{0};

        for (auto && t : text)
        {
            builder.set(prefix_sum);
            prefix_sum += std::ranges::distance(t) + 1;
        }

        text_begin    = sdsl::sd_vector<>(builder);
        text_begin_ss = sdsl::select_support_sd<1>(&text_begin);
        text_begin_rs = sdsl::rank_support_sd<1>(&text_begin);

        // last text in collection needs no delimiter if we have more than one text in the collection
        sdsl::int_vector<8> tmp_text(text_size - (std::ranges::distance(text) > 1));

        constexpr uint8_t delimiter = sigma >= 255 ? 255 : sigma + 1;


        std::ranges::move(text
                          | views::deep{views::to_rank}
                          | views::deep
                          {
                              std::views::transform([] (uint8_t const r)
                              {
                                  if constexpr (sigma >= 255)
                                  {
                                      if (r >= 254)
                                          throw std::out_of_range("The input text cannot be indexed, because"
                                                                  " for full character alphabets the last one/"
                                                                  "two values are reserved (single sequence/"
                                                                  "collection).");
                                  }
                                  return r + 1;
                              })
                          }
                          | views::join(delimiter),
                          std::ranges::begin(tmp_text));


        // we need at least one delimiter
        if (std::ranges::distance(text) == 1)
            tmp_text.back() = delimiter;

        std::ranges::reverse(tmp_text);

        sdsl::construct_im(index, tmp_text, 0);
    }

public:
    //!\brief Indicates whether index is built over a collection.
    static constexpr text_layout text_layout_mode = text_layout_mode_;

    /*!\name Member types
     * \{
     */
    //!\brief The type of the underlying character of the indexed text.
    using alphabet_type = alphabet_t;
    //!\brief Type for representing positions in the indexed text.
    using size_type = typename sdsl_index_type::size_type;
    //!\brief The type of the (unidirectional) cursor.
    using cursor_type = fm_index_cursor<fm_index>;
    //!\}

    template <typename bi_fm_index_t>
    friend class bi_fm_index_cursor;

    template <typename fm_index_t>
    friend class fm_index_cursor;

    template <typename fm_index_t>
    friend class detail::fm_index_cursor_node;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    fm_index() = default;              //!< Defaulted.

    //!\brief When copy constructing, also update internal data structures.
    fm_index(fm_index const & rhs) :
        index{rhs.index}, text_begin{rhs.text_begin}, text_begin_ss{rhs.text_begin_ss}, text_begin_rs{rhs.text_begin_rs}
    {
        text_begin_ss.set_vector(&text_begin);
        text_begin_rs.set_vector(&text_begin);
    }

    //!\brief When move constructing, also update internal data structures.
    fm_index(fm_index && rhs) :
        index{std::move(rhs.index)}, text_begin{std::move(rhs.text_begin)},text_begin_ss{std::move(rhs.text_begin_ss)},
        text_begin_rs{std::move(rhs.text_begin_rs)}
    {
        text_begin_ss.set_vector(&text_begin);
        text_begin_rs.set_vector(&text_begin);
    }

    //!\brief When copy/move assigning, also update internal data structures.
    fm_index & operator=(fm_index rhs)
    {
        index = std::move(rhs.index);
        text_begin = std::move(rhs.text_begin);
        text_begin_ss = std::move(rhs.text_begin_ss);
        text_begin_rs = std::move(rhs.text_begin_rs);

        text_begin_ss.set_vector(&text_begin);
        text_begin_rs.set_vector(&text_begin);

        return *this;
    }

    ~fm_index() = default;             //!< Defaulted.

    /*!\brief Constructor that immediately constructs the index given a range. The range cannot be empty.
     * \tparam text_t The type of range to construct from; must model std::ranges::bidirectional_range.
     * \param[in] text The text to construct from.
     *
     * ### Complexity
     *
     * \if DEV \todo \endif At least linear.
     */
    template <std::ranges::bidirectional_range text_t>
    explicit fm_index(text_t && text)
    {
        construct(std::forward<text_t>(text));
    }
    //!\}

    /*!\brief Returns the length of the indexed text including sentinel characters.
     * \returns Returns the length of the indexed text including sentinel characters.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    size_type size() const noexcept
    {
        return index.size();
    }

    /*!\brief Checks whether the index is empty.
     * \returns `true` if the index is empty, `false` otherwise.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    bool empty() const noexcept
    {
        return size() == 0;
    }

    /*!\brief Compares two indices.
     * \returns `true` if the indices are equal, false otherwise.
     *
     * ### Complexity
     *
     * Linear.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    bool operator==(fm_index const & rhs) const noexcept
    {
        // (void) rhs;
        return (index == rhs.index) && (text_begin == rhs.text_begin);
    }

    /*!\brief Compares two indices.
     * \returns `true` if the indices are unequal, false otherwise.
     *
     * ### Complexity
     *
     * Linear.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    bool operator!=(fm_index const & rhs) const noexcept
    {
        return !(*this == rhs);
    }

    /*!\brief Returns a seqan3::fm_index_cursor on the index that can be used for searching.
     *        \if DEV
     *            Cursor is pointing to the root node of the implicit suffix tree.
     *        \endif
     * \returns Returns a (unidirectional) seqan3::fm_index_cursor on the index.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    cursor_type cursor() const noexcept
    {
        return {*this};
    }

    /*!\cond DEV
     * \brief Serialisation support function.
     * \tparam archive_t Type of `archive`; must satisfy seqan3::cereal_archive.
     * \param archive The archive being serialised from/to.
     *
     * \attention These functions are never called directly, see \ref serialisation for more details.
     */
    template <cereal_archive archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        archive(index);
        archive(text_begin);
        archive(text_begin_ss);
        text_begin_ss.set_vector(&text_begin);
        archive(text_begin_rs);
        text_begin_rs.set_vector(&text_begin);

        auto sigma = alphabet_size<alphabet_t>;
        archive(sigma);
        if (sigma != alphabet_size<alphabet_t>)
        {
            throw std::logic_error{"The fm_index was built over an alphabet of size " + std::to_string(sigma) +
                                   " but it is being read into an fm_index with an alphabet of size " +
                                   std::to_string(alphabet_size<alphabet_t>) + "."};
        }

        bool tmp = text_layout_mode;
        archive(tmp);
        if (tmp != text_layout_mode)
        {
            throw std::logic_error{std::string{"The fm_index was built over a "} +
                                   (tmp ? "text collection" : "single text") +
                                   " but it is being read into an fm_index expecting a " +
                                   (text_layout_mode ? "text collection." : "single text.")};
        }
    }
    //!\endcond

};

/*!\name Template argument type deduction guides
 * \{
 */
//! \brief Deduces the alphabet and dimensions of the text.
template <std::ranges::range text_t>
fm_index(text_t &&) -> fm_index<innermost_value_type_t<text_t>, text_layout{dimension_v<text_t> != 1}>;
//!\}

//!\}
} // namespace seqan3
