// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \brief Provides the unidirectional seqan3::fm_index.
 */

#pragma once

#include <sdsl/suffix_trees.hpp>

#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/range/shortcuts.hpp>
#include <seqan3/range/view/to_rank.hpp>
#include <seqan3/search/fm_index/concept.hpp>
#include <seqan3/search/fm_index/detail/csa_alphabet_strategy.hpp>
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

/*!\brief The default FM Index Configuration.
 *
 * \details
 *
 * ### Running time / Space consumption
 *
 * \f$SAMPLING\_RATE = 16\f$
 * \f$\Sigma\f$: alphabet_size<char_type> where char_type is the seqan3 alphabet type (e.g. dna4 has an alphabet size of 4).
 *
 * For an index over a text collection a delimiter is added inbetween the texts. This causes sigma to increase by 1.
 * \attention For any alphabet, the symbol with rank 255 is not allowed to occur in the text. Addtionally,
 *            rank 254 cannot occur when indexing text collections.
 *
 * \if DEV
 * Rank 255 is generally not allowed. When constructing the index, we increase the rank of every letter by 1, since the
 * SDSL uses 0 as a sentinel. Incrementing 255 by 1 will cause the rank (uint8_t) to overflow to 0.
 * For text collections, 255 is reserved as delimiter between the individual texts. Therefore the letter with rank 254
 * cannot occur in the text.
 *
 * This index will only work for byte alphabets, i.e. alphabets with a size <= 256.
 * When switching to bigger alphabets, the requires clause of the seqan3::fm_index, the seqan3::fm_index_default_traits,
 * the delimiter choice and the transform view for the construction need to be adjusted. Additionally, all occurrences
 * of uint8_t should be double checked to make sure they also apply to bigger alphabets.
 * \endif
 *
 * \f$T_{BACKWARD\_SEARCH}: O(\log \Sigma)\f$
 *
 * \todo Asymptotic space consumption:
 *
 */
struct fm_index_default_traits
{
    //!\brief Type of the underlying SDSL index.
    using sdsl_index_type = sdsl::csa_wt<
        sdsl::wt_blcd<
            sdsl::bit_vector,
            sdsl::rank_support_v<>,
            sdsl::select_support_scan<>,
            sdsl::select_support_scan<0>
        >,
        16,
        10000000,
        sdsl::sa_order_sa_sampling<>,
        sdsl::isa_sampling<>,
        sdsl::plain_byte_alphabet
    >;
};

/*!\brief The SeqAn FM Index.
 * \implements seqan3::FmIndex
 * \tparam text_t The type of the text to be indexed; must model std::ranges::ForwardRange.
 * \tparam fm_index_traits The traits determining the implementation of the underlying SDSL index;
                           must model seqan3::FmIndexTraits.
 * \details
 *
 * The seqan3::fm_index is a fast and space-efficient string index to search strings and collections of strings.
 *
 * ### General information
 *
 * Here is a short example on how to build an index and search a pattern using an cursor. Please note that there is a
 * very powerful search module with a high-level interface \todo seqan3::search that encapsulates the use of cursors.
 *
 * \include test/snippet/search/fm_index.cpp
 *
 * \attention When building an index for a text collection over any alphabet, the symbol with rank 255 is reserved
 *            and may not occur in the text.
 *
 * Here is an example using a collection of strings (e.g. a genome with multiple chromosomes or a protein database):
 *
 * \include test/snippet/search/fm_index_collection.cpp
 *
 * \attention When building an index for a text collection over any alphabet, the symbols with rank 254 and 255
              are reserved and may not be used in the text.
 *
 * ### Choosing an index implementation
 *
 * \todo The underlying implementation of the FM Index (Rank data structure, sampling rates, etc.) can be specified ...
 */
template <std::ranges::RandomAccessRange text_t, FmIndexTraits fm_index_traits = fm_index_default_traits>
//!\cond
    requires Semialphabet<innermost_value_type_t<text_t>> &&
             alphabet_size<innermost_value_type_t<text_t>> <= 256
//!\endcond
class fm_index
{
protected:
    //!\privatesection

    /*!\name Member types
     * \{
     */
    //!\brief The type of the underlying SDSL index.
    using sdsl_index_type = typename fm_index_traits::sdsl_index_type;
    /*!\brief The type of the reduced alphabet type. (The reduced alphabet might be smaller than the original alphabet
     *        in case not all possible characters occur in the indexed text.)
     */
    using sdsl_char_type = typename sdsl_index_type::alphabet_type::char_type;
    //!\brief The type of the alphabet size of the underlying SDSL index.
    using sdsl_sigma_type = typename sdsl_index_type::alphabet_type::sigma_type;
    //!\}

    //!\brief Underlying index from the SDSL.
    sdsl_index_type index;
    //!\brief Pointer to the indexed text.
    text_t const * text = nullptr;

    //!\brief Bitvector storing begin positions for collections.
    sdsl::sd_vector<> text_begin;
    //!\brief Select support for text_begin.
    sdsl::select_support_sd<1> text_begin_ss;
    //!\brief Rank support for text_begin.
    sdsl::rank_support_sd<1> text_begin_rs;

public:
    /*!\name Member types
     * \{
     */
    //!\brief The type of the indexed text.
    using text_type = text_t;
    //!\brief The type of the underlying character of text_type.
    using char_type = innermost_value_type_t<text_t>;
    //!\brief Type for representing positions in the indexed text.
    using size_type = typename sdsl_index_type::size_type;
    //!\brief The type of the (unidirectional) cursor.
    using cursor_type = fm_index_cursor<fm_index<text_t, fm_index_traits>>;
    //!\}

    static_assert(dimension_v<text_t> == 1 || dimension_v<text_t> == 2,
                  "Only texts or collections of texts can be indexed.");

    //!\brief Indicates whether index is built over a collection.
    static bool constexpr is_collection = dimension_v<text_t> == 2;

    template <typename bi_fm_index_t>
    friend class bi_fm_index_cursor;

    template <typename fm_index_t>
    friend class fm_index_cursor;

    template <typename fm_index_t>
    friend class detail::fm_index_cursor_node;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    fm_index() = default;                             //!< Default constructor.
    fm_index(fm_index const &) = default;             //!< Copy constructor.
    fm_index & operator=(fm_index const &) = default; //!< Copy assignment.
    fm_index(fm_index &&) = default;                  //!< Move constructor.
    fm_index & operator=(fm_index &&) = default;      //!< Move assignment.
    ~fm_index() = default;                            //!< Destructor.

    /*!\brief Constructor that immediately constructs the index given a range.
              The range cannot be an rvalue (i.e. a temporary object) and has to be non-empty.
     * \tparam text_t The type of range to construct from; must model std::ranges::RandomAccessRange.
     * \param[in] text The text to construct from.
     *
     * ### Complexity
     *
     * \todo At least linear.
     */
    fm_index(text_t const & text)
    {
        construct(text);
    }

    //!\overload
    fm_index(text_t &&) = delete;

    //!\overload
    fm_index(text_t const &&) = delete;
    //!\}

    /*!\brief Constructs the index given a range.
              The range cannot be an rvalue (i.e. a temporary object) and has to be non-empty.
     * \tparam text_t The type of range to construct from; must model std::ranges::RandomAccessRange.
     * \param[in] text The text to construct from.
     *
     * \details \todo This has to be better implemented with regard to the memory peak due to not matching interfaces
     *                with the SDSL.
     *
     * ### Complexity
     *
     * \todo At least linear.
     *
     * ### Exceptions
     *
     * No guarantees.
     */
    void construct(text_t const & text)
        //!\cond
        requires !is_collection
        //!\endcond
    {
         // text must not be empty
        if (std::ranges::begin(text) == std::ranges::end(text))
            throw std::invalid_argument("The text that is indexed cannot be empty.");

        this->text = &text;
        // TODO:
        // * check what happens in sdsl when constructed twice!
        // * choose between in-memory/external and construction algorithms
        // * sdsl construction currently only works for int_vector, std::string and char *, not ranges in general
        // uint8_t largest_char = 0;
        sdsl::int_vector<8> tmp_text(text.size());

        std::ranges::copy(text
                          | view::to_rank
                          | std::view::transform([] (uint8_t const r)
                          {
                              if constexpr (alphabet_size<char_type> == 256)
                              {
                                  if (r == 255)
                                      throw std::out_of_range("The input text cannot be indexed, because for full"
                                                              "character alphabets the last one/two values are reserved"
                                                              "(single sequence/collection).");
                              }
                              return r + 1;
                          })
                          | std::view::reverse,
                          seqan3::begin(tmp_text)); // reverse and increase rank by one

        sdsl::construct_im(index, tmp_text, 0);

        // TODO: would be nice but doesn't work since it's private and the public member references are const
        // index.m_C.resize(largest_char);
        // index.m_C.shrink_to_fit();
        // index.m_sigma = largest_char;
    }

    //!\overload
    void construct(text_t const & text)
        //!\cond
        requires is_collection
        //!\endcond
    {
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
            text_size += 1 + t.size(); // text size and delimiter (sum will be 1 for empty texts)
        }

        if (all_empty)
            throw std::invalid_argument("A text collection that only contains empty texts cannot be indexed.");

        this->text = &text;

        // bitvector where 1 marks the begin position of a single text from the collection in the concatenated text
        sdsl::bit_vector pos(text_size, 0);
        size_t prefix_sum{0};

        for (auto && t : text)
        {
            pos[prefix_sum] = 1;
            prefix_sum += t.size() + 1;
        }

        text_begin    = sdsl::sd_vector(pos);
        text_begin_ss = sdsl::select_support_sd<1>(&text_begin);
        text_begin_rs = sdsl::rank_support_sd<1>(&text_begin);

        sdsl::int_vector<8> tmp_text(text_size - 1); // last text in collection needs no delimiter

        uint8_t delimiter = alphabet_size<char_type> >= 255 ? 255 : alphabet_size<char_type> + 1;

        std::vector<uint8_t> tmp = text
                                   | view::deep{view::to_rank}
                                   | view::deep
                                   {
                                       std::view::transform([] (uint8_t const r)
                                       {
                                           if constexpr (alphabet_size<char_type> >= 255)
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
                                   | std::view::join(delimiter);

        std::ranges::copy((tmp | std::view::reverse), seqan3::begin(tmp_text));

        //!\todo Replace with this once this does not cause debug builds to exceed max memory on travis
        // std::ranges::copy(text
        //                   | view::deep{view::to_rank}
        //                   | view::deep{std::view::transform([] (uint8_t const r) { return r + 1; })} // increase rank
        //                   | view::deep{std::view::reverse}
        //                   | std::view::reverse
        //                   | std::view::join(delimiter), // join with delimiter
        //                   seqan3::begin(tmp_text));

        sdsl::construct_im(index, tmp_text, 0);
    }

    //!\overload
    void construct(text_t &&) = delete;

    //!\overload
    void construct(text_t const &&) = delete;

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
    cursor_type begin() const noexcept
    {
        return {*this};
    }

    /*!\cond DEV
     * \brief Serialisation support function.
     * \tparam archive_t Type of `archive`; must satisfy seqan3::CerealArchive.
     * \param archive The archive being serialised from/to.
     *
     * \attention These functions are never called directly, see \ref serialisation for more details.
     */
    template <CerealArchive archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        archive(index);
        archive(text_begin);
        archive(text_begin_ss);
        text_begin_ss.set_vector(&text_begin);
        archive(text_begin_rs);
        text_begin_rs.set_vector(&text_begin);
    }
    //!\endcond

};

//!\}

} // namespace seqan3
