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
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \brief Provides the bidirectional seqan3::bi_fm_index.
 */

#pragma once

#include <utility>

#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/io/filesystem.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/fm_index/bi_fm_index_iterator.hpp>
#include <seqan3/std/view/reverse.hpp>

namespace seqan3
{

/*!\addtogroup submodule_fm_index
 * \{
 */

/*!\brief The default Bidirectional FM Index Configuration.
 *
 * \details
 *
 * \todo write me
 *
 */
struct bi_fm_index_default_traits
{
    //!\brief Type of the underlying forward SDSL index.
    using fm_index_traits = fm_index_default_traits;

    //!\brief Type of the underlying reverse SDSL index.
    using rev_fm_index_traits = fm_index_default_traits; // TODO: trait object without sampling.
};

/*!\brief The SeqAn Bidirectional FM Index
 * \implements seqan3::BiFmIndex
 * \tparam text_t The type of the text to be indexed; must model std::ranges::RandomAccessRange.
 * \tparam bi_fm_index_traits The traits determining the implementation of the underlying SDSL indices;
                              must model seqan3::BiFmIndexTraits.
 * \details
 *
 * \todo write me
 */
template <std::ranges::RandomAccessRange text_t, BiFmIndexTraits index_traits_t = bi_fm_index_default_traits>
//!\cond
    requires Alphabet<innermost_value_type_t<text_t>> &&
             std::Same<underlying_rank_t<innermost_value_type_t<text_t>>, uint8_t>
//!\endcond
class bi_fm_index
{
protected:
    //!\privatesection

    //!\brief Pointer to the indexed text.
    text_t const * text = nullptr;

public:

    /*!\name Text types
     * \{
     */
    //!\brief The type of the forward indexed text.
    using text_type = text_t;
    //!\brief The type of the forward indexed text.
    using rev_text_type = decltype(view::reverse(*text));
    //!\}

protected:
    //!\privatesection

    /*!\name Index types
     * \{
     */
    //!\brief The type of the underlying SDSL index for the original text.
    using sdsl_index_type = typename index_traits_t::fm_index_traits::sdsl_index_type;

    //!\brief The type of the underlying SDSL index for the reversed text.
    using rev_sdsl_index_type = typename index_traits_t::rev_fm_index_traits::sdsl_index_type;

    /*!\brief The type of the reduced alphabet type. (The reduced alphabet might be smaller than the original alphabet
     *        in case not all possible characters occur in the indexed text.)
     */
    using sdsl_char_type = typename sdsl_index_type::alphabet_type::char_type;

    //!\brief The type of the underlying FM index for the original text.
    using fm_index_type = fm_index<text_t, typename index_traits_t::fm_index_traits>;

    //!\brief The type of the underlying FM index for the reversed text.
    using rev_fm_index_type = fm_index<rev_text_type, typename index_traits_t::rev_fm_index_traits>;
    //!\}

    //!\brief Access to a reversed view of the text. Needed when a unidirectional iterator on the reversed text is
    //        constructed from the bidirectional index.
    rev_text_type rev_text;

    //!\brief Underlying FM index for the original text.
    fm_index_type fwd_fm;

    //!\brief Underlying FM index for the reversed text.
    rev_fm_index_type rev_fm;

public:

    /*!\name Text types
     * \{
     */
    //!\brief The type of the underlying character of text_type.
    using char_type = innermost_value_type_t<text_t>;
    //!\brief Type for representing positions in the indexed text.
    using size_type = typename sdsl_index_type::size_type;
    //!\}

    //!\brief The index traits object.
    using index_traits = index_traits_t;

    /*!\name Iterator types
     * \{
     */
    //!\brief The type of the bidirectional iterator.
    using iterator_type = bi_fm_index_iterator<bi_fm_index<text_t, index_traits_t>>;
    //!\brief The type of the unidirectional iterator on the original text.
    using fwd_iterator_type = fm_index_iterator<fm_index_type>;
    //!\brief The type of the unidirectional iterator on the reversed text.
    using rev_iterator_type = fm_index_iterator<rev_fm_index_type>;
    //!\}

    template <typename bi_fm_index_t>
    friend class bi_fm_index_iterator;

    template <typename fm_index_t>
    friend class fm_index_iterator;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    bi_fm_index() = default;
    bi_fm_index(bi_fm_index const &) = default;
    bi_fm_index & operator=(bi_fm_index const &) = default;
    bi_fm_index(bi_fm_index &&) = default;
    bi_fm_index & operator=(bi_fm_index &&) = default;
    ~bi_fm_index() = default;

    /*!\brief Constructor that immediately constructs the index given a range.
     *        The range cannot be an rvalue (i.e. a temporary object) and has to be non-empty.
     * \tparam text_t The type of range to construct from; must model std::ranges::RandomAccessRange.
     * \param[in] text The text to construct from.
     *
     * ### Complexity
     *
     * \todo At least linear.
     */
    bi_fm_index(text_t const & text)
    {
        construct(text);
    }

    //!\overload
    bi_fm_index(text_t &&) = delete;

    //!\overload
    bi_fm_index(text_t const &&) = delete;
    //!\}

    /*!\brief Constructs the index given a range.
     *        The range cannot be an rvalue (i.e. a temporary object) and has to be non-empty.
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
     * No guarantee. \todo Ensure strong exception guarantee.
     */
    void construct(text_t const & text)
    {
         // text must not be empty
        if (text.begin() == text.end())
            throw std::invalid_argument("The text that is indexed cannot be empty.");

        this->text = &text;
        rev_text = view::reverse(text);
        fwd_fm.construct(text);
        rev_fm.construct(rev_text);

        // does not work yet. segmentation fault in bi_fm_index_iterator snippet
        // bi_fm_index tmp;
        // tmp.text = &text;
        // tmp.rev_text = view::reverse(*tmp.text);
        // tmp.fwd_fm.construct(*tmp.text);
        // tmp.rev_fm.construct(tmp.rev_text);
        // std::swap(*this, tmp);
        // this->text = &text;
        // rev_text = view::reverse(text);
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

    // operator== not implemented by sdsl indices yet
    // bool operator==(fm_index const & rhs) const noexcept
    // {
    //     return std::tie(fwd_fm, rev_fm) == std::tie(rhs.fwd_fm, rhs.rev_fm);
    // }

    // operator== not implemented by sdsl indices yet
    // bool operator!=(fm_index const & rhs) const noexcept
    // {
    //     return !(*this == rhs);
    // }

    /*!\brief Returns a seqan3::bi_fm_index_iterator on the index that can be used for searching.
     *        \cond DEV
     *            Iterator is pointing to the root node of the implicit affix tree.
     *        \endcond
     * \returns Returns a bidirectional seqan3::bi_fm_index_iterator on the index.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    iterator_type begin() const noexcept
    {
        return {*this};
    }

    /*!\brief Returns a unidirectional seqan3::fm_index_iterator on the original text of the bidirectional index that
     *        can be used for searching.
     * \returns Returns a unidirectional seqan3::fm_index_iterator on the index of the original text.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    fwd_iterator_type fwd_begin() const noexcept
    {
       return {fwd_fm};
    }

    /*!\brief Returns a unidirectional seqan3::fm_index_iterator on the reversed text of the bidirectional index that
     *        can be used for searching. Note that because of the text being reversed, extend_right() resp. cycle_back()
     *        correspond to extend_left() resp. cycle_front() on the bidirectional index iterator.
     * \returns Returns a unidirectional seqan3::fm_index_iterator on the index of the reversed text.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    rev_iterator_type rev_begin() const noexcept
    {
       return {rev_fm};
    }

    /*!\brief Loads the index from disk. Temporary function until cereal is supported.
     *        \todo cereal
     * \returns `true` if the index was successfully loaded from disk.
     *
     * ### Complexity
     *
     * Linear.
     *
     * ### Exceptions
     *
     * Strong exception guarantee.
     */
    bool load(filesystem::path const & path)
    {
        filesystem::path path_fwd{path};
        filesystem::path path_rev{path};
        path_fwd += filesystem::path{".fwd"};
        path_rev += filesystem::path{".rev"};
        return fwd_fm.load(path_fwd) && rev_fm.load(path_rev);
    }

    /*!\brief Stores the index to disk. Temporary function until cereal is supported.
     *        \todo cereal
     * \returns `true` if the index was successfully stored to disk.
     *
     * ### Complexity
     *
     * Linear.
     *
     * ### Exceptions
     *
     * Strong exception guarantee.
     */
    bool store(filesystem::path const & path) const
    {
        filesystem::path path_fwd{path};
        filesystem::path path_rev{path};
        path_fwd += filesystem::path{".fwd"};
        path_rev += filesystem::path{".rev"};
        return fwd_fm.store(path_fwd) && rev_fm.store(path_rev);
    }

};

//!\}

} // namespace seqan3
