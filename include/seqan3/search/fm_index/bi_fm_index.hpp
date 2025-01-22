// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \brief Provides the bidirectional seqan3::bi_fm_index.
 */

#pragma once

#include <filesystem>
#include <ranges>
#include <utility>

#include <seqan3/core/range/type_traits.hpp>
#include <seqan3/search/fm_index/bi_fm_index_cursor.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>

namespace seqan3
{

/*!\brief The SeqAn Bidirectional FM Index
 * \ingroup search_fm_index
 * \tparam alphabet_t        The alphabet type; must model seqan3::semialphabet.
 * \tparam text_layout_mode_ Indicates whether this index works on a text collection or a single text.
 *                           See seqan3::text_layout.
 * \tparam sdsl_index_type_  The type of the underlying SDSL index, must model seqan3::detail::sdsl_index.
 * \implements seqan3::cerealisable
 * \details
 *
 * The seqan3::bi_fm_index is a fast and space-efficient bidirectional string index to search strings and
 * collections of strings.
 * In general, we recommend to favour the seqan3::bi_fm_index over the unidirectional seqan3::fm_index if you want to
 * allow multiple errors when searching.
 *
 * ### General information
 *
 * Here is a short example on how to build an index and search a pattern using an cursor. Please note that there is a
 * very powerful search module with a high-level interface seqan3::search that encapsulates the use of cursors.
 *
 * \include test/snippet/search/bi_fm_index.cpp
 *
 * \attention When building an index for a **single text** over any alphabet, the symbol with rank 255 is reserved
 *            and may not occur in the text.
 *
 * Here is an example using a collection of strings (e.g. a genome with multiple chromosomes or a protein database):
 *
 * \include test/snippet/search/bi_fm_index_collection.cpp
 *
 * \attention When building an index for a **text collection** over any alphabet, the symbols with rank 254 and 255
 *            are reserved and may not be used in the text.
 */
template <semialphabet alphabet_t,
          text_layout text_layout_mode_,
          detail::sdsl_index sdsl_index_type_ = default_sdsl_index_type>
class bi_fm_index
{
private:
    /*!\name Index types
     * \{
     */
    //!\brief The type of the underlying SDSL index for the original text.
    using sdsl_index_type = sdsl_index_type_;

    //!\brief The type of the underlying SDSL index for the reversed text.
    using rev_sdsl_index_type = sdsl::csa_wt<sdsl_wt_index_type::wavelet_tree_type, // Wavelet tree type
                                             10'000'000,                            // Sampling rate of the suffix array
                                             10'000'000,                   // Sampling rate of the inverse suffix array
                                             sdsl::sa_order_sa_sampling<>, // Text or SA based sampling for SA
                                             sdsl::isa_sampling<>,         // Text or ISA based sampling for ISA
                                             sdsl_wt_index_type::alphabet_type>; // How to represent the alphabet

    /*!\brief The type of the reduced alphabet type. (The reduced alphabet might be smaller than the original alphabet
     *        in case not all possible characters occur in the indexed text.)
     */
    using sdsl_char_type = typename sdsl_index_type::alphabet_type::char_type;

    //!\brief The type of the alphabet size of the underlying SDSL index.
    using sdsl_sigma_type = typename sdsl_index_type::alphabet_type::sigma_type;

    //!\brief The type of the underlying FM index for the original text.
    using fm_index_type = fm_index<alphabet_t, text_layout_mode_, sdsl_index_type>;

    //!\brief The type of the underlying FM index for the reversed text.
    using rev_fm_index_type = detail::reverse_fm_index<alphabet_t, text_layout_mode_, rev_sdsl_index_type>;
    //!\}

    //!\brief Underlying FM index for the original text.
    fm_index_type fwd_fm;

    //!\brief Underlying FM index for the reversed text.
    rev_fm_index_type rev_fm;

    /*!\brief Constructs the index given a range.
     *        The range cannot be an rvalue (i.e. a temporary object) and has to be non-empty.
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
    void construct(text_t && text)
    {
        detail::fm_index_validator::validate<alphabet_t, text_layout_mode_>(text);

        fwd_fm = fm_index_type{text};
        rev_fm = rev_fm_index_type{text};
    }

public:
    //!\brief Indicates whether index is built over a collection.
    static constexpr text_layout text_layout_mode = text_layout_mode_;

    /*!\name Text types
     * \{
     */
    //!\brief The type of the underlying character of the indexed text.
    using alphabet_type = typename fm_index_type::alphabet_type;
    //!\brief Type for representing positions in the indexed text.
    using size_type = typename sdsl_index_type::size_type;
    //!\}

    /*!\name Cursor types
     * \{
     */
    //!\brief The type of the bidirectional cursor.
    using cursor_type = bi_fm_index_cursor<bi_fm_index>;
    //!\brief The type of the unidirectional cursor on the original text.
    using fwd_cursor_type = fm_index_cursor<fm_index_type>;

    //!\}

    template <typename bi_fm_index_t>
    friend class bi_fm_index_cursor;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    bi_fm_index() = default;                                //!< Defaulted.
    bi_fm_index(bi_fm_index const &) = default;             //!< Defaulted.
    bi_fm_index & operator=(bi_fm_index const &) = default; //!< Defaulted.
    bi_fm_index(bi_fm_index &&) = default;                  //!< Defaulted.
    bi_fm_index & operator=(bi_fm_index &&) = default;      //!< Defaulted.
    ~bi_fm_index() = default;                               //!< Defaulted.

    /*!\brief Constructor that immediately constructs the index given a range. The range cannot be empty.
     * \tparam text_t The type of range to construct from; must model std::ranges::bidirectional_range.
     * \param[in] text The text to construct from.
     *
     * ### Complexity
     *
     * \if DEV \todo \endif At least linear.
     */
    template <std::ranges::range text_t>
    bi_fm_index(text_t && text)
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
    cursor_type cursor() const noexcept
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
    fwd_cursor_type fwd_cursor() const noexcept
    {
        return {fwd_fm};
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
        archive(fwd_fm);
        archive(rev_fm);
    }
    //!\endcond
};

/*!\name Template argument type deduction guides
 * \{
 */
//!\brief Deduces the dimensions of the text.
template <std::ranges::range text_t>
bi_fm_index(text_t &&) -> bi_fm_index<range_innermost_value_t<text_t>, text_layout{range_dimension_v<text_t> != 1}>;
//!\}

} // namespace seqan3
