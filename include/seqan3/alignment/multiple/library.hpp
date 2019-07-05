// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 * \brief Provides seqan3::msa_library.
 */

#pragma once

#include <cstddef>
#include <map>
#include <optional>
#include <ostream>
#include <stdexcept>

#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/range/container/concept.hpp>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

/*!\brief Library for multiple sequence alignments.
 * \tparam score_type The type of the score values, must model std::Arithmetic.
 *
 * \details
 * The library stores for each sequence pair a map of positions with associated scores.
 * This class provides accessors, modifiers and an iterator to the underlying data.
 */
template <Arithmetic score_type>
class msa_library
{
public:
    //!\brief The type of an index pair for sequences or positions.
    using coord_type = std::pair<size_t, size_t>;
    //!\brief The type for the map of scores for each position in a pairwise alignment.
    using map_pos_score_type = std::map<coord_type, score_type>;

private:
    //!\brief The iterator type for a sequence pair map.
    using sequence_pair_iterator = typename std::map<coord_type, map_pos_score_type>::iterator;
    //!\brief The iterator type for a position pair map.
    using position_pair_iterator = typename map_pos_score_type::iterator;

    //!\brief A map where each sequence pair is assigned a map of position pairs and scores.
    std::map<coord_type, map_pos_score_type> data;

    /*!\brief Swap the order of both pairs if the first sequence index is greater than the second (to avoid duplicates).
     * \param[in,out] seq A pair of sequence indices.
     * \param[in,out] pos A pair of position indices.
     */
    void check_seq_order(coord_type & seq, coord_type & pos) const noexcept
    {
        if (seq.first > seq.second)
        {
            seq = {seq.second, seq.first};
            pos = {pos.second, pos.first};
        }
    }

    /*!\brief Stream operator for the library.
     * \tparam score_t The score type.
     * \param stream   The stream where the library should be printed to.
     * \return         The modified stream object.
     */
    template <typename score_t>
    friend std::ostream & operator<<(std::ostream & stream, msa_library<score_t> const &);

    /*!\brief The iterator for seqan3::msa_library.
     *        The elements are accessed in the order of sequence and position indices.
     */
    class msa_library_iterator
    {
    public:
        //!\brief The reference type. As sequence and position pairs are constant, only the score is modifiable.
        using reference = std::tuple<coord_type const &, coord_type const &, score_type &>;
        //!\brief The value type is a std::tuple of sequence pair, position pair and the score.
        using value_type = std::tuple<coord_type const, coord_type const, score_type>;
        //!\brief The difference type of iterators.
        using difference_type = std::ptrdiff_t;

    private:
        //!\brief Pointer to the underlying msa library.
        typename std::add_pointer_t<msa_library> host{nullptr};
        //!\brief Internal iterator for the sequence pair map.
        sequence_pair_iterator seq_it{};
        //!\brief Internal iterator for the position pair maps.
        position_pair_iterator pos_it{};

    public:
        /*!\name Constructors/Destructors
         * \{
         */
        constexpr msa_library_iterator()                                             = default; //!< Defaulted.
        constexpr msa_library_iterator(msa_library_iterator const &)                 = default; //!< Defaulted.
        constexpr msa_library_iterator (msa_library_iterator &&) noexcept            = default; //!< Defaulted.
        constexpr msa_library_iterator & operator=(msa_library_iterator const &)     = default; //!< Defaulted.
        constexpr msa_library_iterator & operator=(msa_library_iterator &&) noexcept = default; //!< Defaulted.
        ~msa_library_iterator()                                                      = default; //!< Defaulted.

        /*!\brief Construct from seqan3::msa_library and set the inner iterator to the specified location.
         * \param lib A pointer to the underlying seqan3::msa_library.
         * \param s_it A sequence pair iterator to specify the location in the outer map.
         * \param p_it A position pair iterator to specify the location in the inner map.
         */
        explicit constexpr msa_library_iterator(msa_library * lib,
                                                sequence_pair_iterator const & s_it,
                                                position_pair_iterator const & p_it) noexcept :
            host(lib),
            seq_it(s_it),
            pos_it(p_it)
        {}

        /*!\brief Construct from seqan3::msa_library and set members to the begin state.
         * \param lib A pointer to the underlying seqan3::msa_library.
         */
        explicit constexpr msa_library_iterator(msa_library * lib) noexcept :
            msa_library_iterator(lib, lib->data.begin(), {})
        {
            if (seq_it != host->data.end())
                pos_it = seq_it->second.begin();
        }
        //!\}

        /*!\name Arithmetic operators
         * \{
        */
        //!\brief Pre-increment, returns updated iterator.
        constexpr msa_library_iterator & operator++() noexcept
        {
            // We have to test whether seq_it is at end, because we may not be able to dereference it.
            if (seq_it != host->data.end() && ++pos_it == seq_it->second.end() && ++seq_it != host->data.end())
                pos_it = seq_it->second.begin();
            return *this;
        }

        //!\brief Pre-decrement, returns updated iterator.
        constexpr msa_library_iterator & operator--() noexcept
        {
            // We have to test whether seq_it is invalid, because we may not be able to dereference it.
            if ((seq_it == host->data.end() || pos_it == seq_it->second.begin()) && seq_it != host->data.begin())
                pos_it = (--seq_it)->second.end();

            --pos_it;
            return *this;
        }

        //!\brief Post-increment, returns previous iterator state.
        constexpr msa_library_iterator operator++(int) noexcept
        {
            msa_library_iterator cpy{*this};
            ++(*this);
            return cpy;
        }

        //!\brief Post-decrement, returns previous iterator state.
        constexpr msa_library_iterator operator--(int) noexcept
        {
            msa_library_iterator cpy{*this};
            --(*this);
            return cpy;
        }
        //!\}

        /*!\name Access to the data members.
         * \{
         */
        //!\brief Dereference operator returns a std::tuple of sequence pair, position pair and score reference.
        constexpr reference operator*() const noexcept
        {
            // throw std::out_of_range{"Trying to access the element behind the last one in msa_library."};
            assert(seq_it != host->data.end());
            return std::tie(seq_it->first, pos_it->first, pos_it->second);
        }

        //!\brief Returns the score of the current iterator.
        constexpr score_type & score() const noexcept
        {
            assert(seq_it != host->data.end());
            return pos_it->second;
        }

        //!\brief Returns the pair of sequence indices of the current iterator.
        constexpr coord_type const & seq_pair() const noexcept
        {
            assert(seq_it != host->data.end());
            return seq_it->first;
        }

        //!\brief Returns the pair of position indices of the current iterator.
        constexpr coord_type const & pos_pair() const noexcept
        {
            assert(seq_it != host->data.end());
            return pos_it->first;
        }
        //!\}

        /*!\name Comparison operators
         * \brief Compares iterators.
         * \{
         */

        //!\brief Checks whether `*this` is equal to `rhs`.
        constexpr friend bool operator==(msa_library_iterator const & lhs, msa_library_iterator const & rhs)
        {
            return lhs.seq_it == rhs.seq_it && lhs.pos_it == rhs.pos_it;
        }

        //!\brief Checks whether `*this` is not equal to `rhs`.
        constexpr friend bool operator!=(msa_library_iterator const & lhs, msa_library_iterator const & rhs)
        {
            return !(lhs == rhs);
        }

        //!\brief Checks whether `*this` is less than `rhs`.
        constexpr friend bool operator<(msa_library_iterator const & lhs, msa_library_iterator const & rhs)
        {
            return lhs.seq_it < rhs.seq_it || (lhs.seq_it == rhs.seq_it && lhs.pos_it < rhs.pos_it);
        }

        //!\brief Checks whether `*this` is less than or equal to `rhs`.
        constexpr friend bool operator<=(msa_library_iterator const & lhs, msa_library_iterator const & rhs)
        {
            return lhs < rhs || lhs == rhs;
        }

        //!\brief Checks whether `*this` is greater than `rhs`.
        constexpr friend bool operator>(msa_library_iterator const & lhs, msa_library_iterator const & rhs)
        {
            return !(lhs <= rhs);
        }

        //!\brief Checks whether `*this` is greater than or equal to `rhs`.
        constexpr friend bool operator>=(msa_library_iterator const & lhs, msa_library_iterator const & rhs)
        {
            return !(lhs < rhs);
        }
        //!\}
    }; // class msa_library_iterator

public:
    //!\brief The bidirectional iterator that allows to visit all entries of the library.
    using iterator = msa_library_iterator;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr msa_library()                                      = default; //!< Defaulted.
    constexpr msa_library(msa_library const & other)             = default; //!< Defaulted.
    constexpr msa_library(msa_library &&) noexcept               = default; //!< Defaulted.
    constexpr msa_library & operator=(msa_library const & other) = default; //!< Defaulted.
    constexpr msa_library & operator=(msa_library &&) noexcept   = default; //!< Defaulted.
    ~msa_library()                                               = default; //!< Defaulted.
    //!\}

    /*!\brief Insert a new entry into the library.
     * \param seq A pair of sequence indeices.
     * \param pos A pair of position indices.
     * \param score The assigned score for the specified positions in the specified sequences.
     * \return True if the insertion took place, False if a score has already been assigned.
     */
    bool insert(coord_type seq, coord_type pos, score_type score)
    {
        check_seq_order(seq, pos);

        // try to insert a new sequence pair
        auto [iter, inserted] = data.try_emplace(seq, map_pos_score_type{{pos, score}});

        if (inserted)    // we have added a new sequence pair with the entry
            return true;
        else             // append the score to the existing sequences and return whether it was successful
            return iter->second.try_emplace(pos, score).second;
    }

    /*!\brief Add a score to a possibly existing entry in the library.
     * \param seq A pair of sequence indeices.
     * \param pos A pair of position indices.
     * \param score The score that is added to the specified position in the specified sequence pair.
     *              If this location does not exist it will be created with value `score`.
     */
    void add(coord_type seq, coord_type pos, score_type score)
    {
        check_seq_order(seq, pos);

        // try to insert a new sequence pair
        auto [iter, inserted] = data.try_emplace(seq, map_pos_score_type{{pos, score}});

        if (!inserted)   // append to an existing sequence pair
            iter->second[pos] += score;
    }

    /*!
     * \brief Retrieve the iterator that contains the score of pairing two positions in the specified sequences.
     * \param seq_pos A quadruple consisting of (seq1, seq2, pos1, pos2).
     * \return An iterator to the score if the specified entry exists and end() otherwise.
     */
    iterator operator[](std::tuple<size_t, size_t, size_t, size_t> seq_pos) noexcept
    {
        coord_type seq = {std::get<0>(seq_pos), std::get<1>(seq_pos)};
        coord_type pos = {std::get<2>(seq_pos), std::get<3>(seq_pos)};
        check_seq_order(seq, pos);
        auto seq_it = data.find(seq);
        if (seq_it != data.end())
        {
            auto pos_it = seq_it->second.find(pos);
            if (pos_it != seq_it->second.end())
                return iterator(this, seq_it, pos_it);
        }
        return end();
    }

    /*!\brief Retrieve a map of all scores in the two specified sequences.
     * \param seq A pair of sequence indices. The smaller index has to be in the first position of the pair.
     * \return A std::optional that contains a score map if the specified sequences exist and False otherwise.
     */
    std::optional<map_pos_score_type> operator[](std::pair<size_t, size_t> seq) const noexcept
    {
        assert(seq.first < seq.second); // cannot swap the order without rewriting the resulting map

        auto seq_it = data.find(seq);
        if (seq_it != data.end())
            return seq_it->second;
        else
            return {};
    }

    /*!\name Iterators
     * \{
     */
    /*!\brief Returns an iterator to the first element of the library.
     * \returns Iterator to the first element.
     *
     * If the container is empty, the returned iterator will be equal to end().
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    constexpr iterator begin() noexcept
    {
        return iterator(this);
    }

    /*!\brief Returns an iterator to the element following the last element of the library.
     * \returns Iterator behind the last element.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    constexpr iterator end() noexcept
    {
        typename map_pos_score_type::iterator pos_it{};
        if (!data.empty())
            pos_it = (--data.end())->second.end();

        return iterator(this, data.end(), pos_it);
    }
};

/*!\brief Stream operator for the msa_library.
 * \tparam score_type The score type.
 * \param stream      The stream where the library should be printed to.
 * \param lib         The library object.
 * \return            The modified stream object.
 *
 * \details
 * Usage:
 * ~~~{.cpp}
 * stream << lib;
 * ~~~
 */
template <typename score_type>
std::ostream & operator<<(std::ostream & stream, msa_library<score_type> const & lib)
{
    for (auto const & [seq, val] : lib.data)
    {
        stream << "# " << seq.first << ' ' << seq.second << std::endl;
        for (auto const & [pos, score] : val)
            stream << pos.first << ' ' << pos.second << ' ' << score << std::endl;
    }
    return stream;
}

/*!\brief Stream operator for msa_library that produces the T-Coffee library format.
 * \tparam score_type The score type.
 * \tparam id_t       The type of the container that contains the ids.
 * \tparam seq_t      The type of the container that contains the sequences.
 * \param stream      The stream where the library, ids and sequences should be printed to.
 * \param obj         A triple consisting of the library object and two containers of equal length,
 *                    where the id and sequence information is stored.
 *
 * \details
 * Usage:
 * ~~~{.cpp}
 * stream << std::make_tuple(lib, ids, seqs);
 * ~~~
 */
template <typename score_type, Container id_t, Container seq_t>
    requires std::ranges::SizedRange<typename seq_t::value_type>
std::ostream & operator<<(std::ostream & stream, std::tuple<msa_library<score_type>, id_t, seq_t> const & obj)
{
    auto const & [lib, ids, seqs] = obj;
    stream << "! T-COFFEE_LIB_FORMAT_01" << std::endl << ranges::size(seqs) << std::endl;

    for (auto && [id_elem, seq_elem] : std::view::zip(ids, seqs))
        stream << id_elem << ' ' << ranges::size(seq_elem) << ' ' << seq_elem << std::endl;

    stream << lib << "! SEQ_1_TO_N" << std::endl;
    return stream;
}

} // namespace seqan3::detail
