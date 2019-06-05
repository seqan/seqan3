// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \brief Provides the bidirectional seqan3::bi_fm_index.
 */

#pragma once

#include <utility>

#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/range/view/persist.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/fm_index/bi_fm_index_cursor.hpp>
#include <seqan3/std/ranges>

namespace seqan3
{

/*!\addtogroup submodule_fm_index
 * \{
 */

/*!\brief The SeqAn Bidirectional FM Index
 * \implements seqan3::BiFmIndex
 * \tparam is_collection    Indicates whether this index works on a text collection (`true`) or a single text (`false`).
 * \tparam sdsl_index_type_ The type of the underlying SDSL index, must model seqan3::SdslIndex.
 * \details
 *
 * \todo write me
 */
template <bool is_collection = false, detail::SdslIndex sdsl_index_type_ = default_sdsl_index_type>
class bi_fm_index
{
protected:
    //!\brief The alphabet size of the text.
    size_t sigma{0};
    //!\brief Indicates whether index is built over a collection.
    static constexpr bool is_collection_{is_collection};

    /*!\name Index types
     * \{
     */
    //!\brief The type of the underlying SDSL index for the original text.
    using sdsl_index_type = sdsl_index_type_;

    //!\brief The type of the underlying SDSL index for the reversed text.
    using rev_sdsl_index_type = sdsl_index_type_;

    /*!\brief The type of the reduced alphabet type. (The reduced alphabet might be smaller than the original alphabet
     *        in case not all possible characters occur in the indexed text.)
     */
    using sdsl_char_type = typename sdsl_index_type::alphabet_type::char_type;

    //!\brief The type of the alphabet size of the underlying SDSL index.
    using sdsl_sigma_type = typename sdsl_index_type::alphabet_type::sigma_type;

    //!\brief The type of the underlying FM index for the original text.
    using fm_index_type = fm_index<is_collection_, sdsl_index_type>;

    //!\brief The type of the underlying FM index for the reversed text. \todo Overwrite sampling behaviour.
    using rev_fm_index_type = fm_index<is_collection_, sdsl_index_type>;
    //!\}

    //!\brief Underlying FM index for the original text.
    fm_index_type fwd_fm;

    //!\brief Underlying FM index for the reversed text.
    rev_fm_index_type rev_fm;

public:
    /*!\name Text types
     * \{
     */
    //!\brief Type for representing positions in the indexed text.
    using size_type = typename sdsl_index_type::size_type;
    //!\}

    /*!\name Cursor types
     * \{
     */
    //!\brief The type of the bidirectional cursor.
    using cursor_type = bi_fm_index_cursor<bi_fm_index<is_collection_, sdsl_index_type>>;
    //!\brief The type of the unidirectional cursor on the original text.
    using fwd_cursor_type = fm_index_cursor<fm_index_type>;
    //!\brief The type of the unidirectional cursor on the reversed text.
    using rev_cursor_type = fm_index_cursor<rev_fm_index_type>;
    //!\}

    template <typename fm_index_t>
    friend class fm_index_cursor;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    bi_fm_index() = default;                                //!< Default constructor.
    bi_fm_index(bi_fm_index const &) = default;             //!< Copy constructor.
    bi_fm_index & operator=(bi_fm_index const &) = default; //!< Copy assignment.
    bi_fm_index(bi_fm_index &&) = default;                  //!< Move constructor.
    bi_fm_index & operator=(bi_fm_index &&) = default;      //!< Move assignment.
    ~bi_fm_index() = default;                               //!< Destructor.

    /*!\brief Constructor that immediately constructs the index given a range. The range cannot be empty.
     * \tparam text_t The type of range to construct from; must model std::ranges::BidirectionalRange.
     * \param[in] text The text to construct from; must model std::ranges::BidirectionalRange.
     *
     * ### Complexity
     *
     * \todo At least linear.
     */
    template <std::ranges::Range text_t>
    bi_fm_index(text_t && text)
    {
        construct(std::forward<decltype(text)>(text));
    }
    //!\}

    /*!\brief Constructs the index given a range.
     *        The range cannot be an rvalue (i.e. a temporary object) and has to be non-empty.
     * \tparam text_t The type of range to construct from; must model std::ranges::BidirectionalRange.
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
     * No guarantee. \todo Ensure strong exception guarantee.
     */
    template <std::ranges::Range text_t>
        //!\cond
        requires !is_collection_
        //!\endcond
    void construct(text_t && text)
    {
        static_assert(std::ranges::BidirectionalRange<text_t>, "The text must be a BidirectionalRange.");
        static_assert(alphabet_size<innermost_value_type_t<text_t>> <= 256, "The alphabet is too big.");
        static_assert(dimension_v<text_t> == 1, "The input cannot be a text collection.");

        // text must not be empty
        if (std::ranges::begin(text) == std::ranges::end(text))
            throw std::invalid_argument("The text that is indexed cannot be empty.");

        auto rev_text = std::view::reverse(text);
        fwd_fm.construct(text);
        rev_fm.construct(rev_text);

        sigma = fwd_fm.sigma;
    }

    //!\overload
    template <std::ranges::Range text_t>
        //!\cond
        requires is_collection_
        //!\endcond
    void construct(text_t && text)
    {
        static_assert(std::ranges::BidirectionalRange<text_t>, "The text must be a BidirectionalRange.");
        static_assert(std::ranges::BidirectionalRange<reference_t<text_t>>,
                      "The elements of the text collection must be a BidirectionalRange.");
        static_assert(alphabet_size<innermost_value_type_t<text_t>> <= 256, "The alphabet is too big.");
        static_assert(dimension_v<text_t> == 2, "The input must be a text collection.");

        // text must not be empty
        if (std::ranges::begin(text) == std::ranges::end(text))
            throw std::invalid_argument("The text that is indexed cannot be empty.");

        auto rev_text = text | view::deep{std::view::reverse} | std::view::reverse;
        fwd_fm.construct(text);
        rev_fm.construct(rev_text);

        sigma = fwd_fm.sigma;
    }

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
        return fwd_fm.size();
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
    bool operator==(bi_fm_index const & rhs) const noexcept
    {
        return std::tie(fwd_fm, rev_fm) == std::tie(rhs.fwd_fm, rhs.rev_fm);
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
    bool operator!=(bi_fm_index const & rhs) const noexcept
    {
        return !(*this == rhs);
    }

    /*!\brief Returns a seqan3::bi_fm_index_cursor on the index that can be used for searching.
     *        \if DEV
     *            Cursor is pointing to the root node of the implicit affix tree.
     *        \endif
     * \returns Returns a bidirectional seqan3::bi_fm_index_cursor on the index.
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

    /*!\brief Returns a unidirectional seqan3::fm_index_cursor on the original text of the bidirectional index that
     *        can be used for searching.
     * \returns Returns a unidirectional seqan3::fm_index_cursor on the index of the original text.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    fwd_cursor_type fwd_begin() const noexcept
    {
       return {fwd_fm};
    }

    /*!\brief Returns a unidirectional seqan3::fm_index_cursor on the reversed text of the bidirectional index that
     *        can be used for searching. Note that because of the text being reversed, extend_right() resp. cycle_back()
     *        correspond to extend_left() resp. cycle_front() on the bidirectional index cursor.
     * \attention For text collections the text IDs are also reversed.
     * \returns Returns a unidirectional seqan3::fm_index_cursor on the index of the reversed text.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    rev_cursor_type rev_begin() const noexcept
    {
       return {rev_fm};
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
        archive(fwd_fm);
        archive(rev_fm);
    }
    //!\endcond
};

/*!\name Template argument type deduction guides
 * \{
 */
//! \brief Deduces the dimensions of the text.
template <std::ranges::Range text_t>
bi_fm_index(text_t &&) -> bi_fm_index<dimension_v<text_t> != 1>;
//!\}

//!\}

} // namespace seqan3
