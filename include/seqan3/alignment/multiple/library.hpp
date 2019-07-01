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

#include <map>
#include <optional>
#include <ostream>

#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/range/container/concept.hpp>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

template <Arithmetic score_type>
class msa_library
{
private:
    //!\brief The type of an index pair for sequences or positions.
    using coord_type = std::pair<size_t, size_t>;
    //!\brief The type for the scores of pairwise alignment.
    using pos_score_type = std::map<coord_type, score_type>;

    /*!\brief Stream operator for the library.
     * \tparam score_t The score type.
     * \param stream   The stream where the library should be printed to.
     * \return         The modified stream object.
     */
    template <typename score_t>
    friend std::ostream & operator<<(std::ostream & stream, msa_library<score_t> const &);

    //!\brief A map where each sequence pair is assigned a map of position pairs and scores.
    std::map<coord_type, pos_score_type> data;

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

public:
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
        auto [iter, inserted] = data.try_emplace(seq, pos_score_type{{pos, score}});

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
        auto [iter, inserted] = data.try_emplace(seq, pos_score_type{{pos, score}});

        if (!inserted)   // append to an existing sequence pair
            iter->second[pos] += score;
    }

    /*!
     * \brief Retrieve the score of pairing two positions in the specified sequences.
     * \param seq_pos A quadruple consisting of (seq1, seq2, pos1, pos2).
     * \return A std::optional that contains the score if the specified entry exists and False otherwise.
     */
    std::optional<score_type> operator[](std::tuple<size_t, size_t, size_t, size_t> seq_pos) const noexcept
    {
        coord_type seq = {std::get<0>(seq_pos), std::get<1>(seq_pos)};
        coord_type pos = {std::get<2>(seq_pos), std::get<3>(seq_pos)};
        check_seq_order(seq, pos);
        auto seq_it = data.find(seq);
        if (seq_it != data.end())
        {
            auto pos_it = seq_it->second.find(pos);
            if (pos_it != seq_it->second.end())
                return pos_it->second;
        }
        return {};
    }

    /*!\brief Retrieve a map of all scores in the two specified sequences.
     * \param seq A pair of sequence indices. The smaller index has to be in the first position of the pair.
     * \return A std::optional that contains a score map if the specified sequences exist and False otherwise.
     */
    std::optional<pos_score_type> operator[](std::pair<size_t, size_t> seq) const noexcept
    {
        assert(seq.first < seq.second); // cannot swap the order without rewriting the resulting map

        auto seq_it = data.find(seq);
        if (seq_it != data.end())
            return seq_it->second;
        else
            return {};
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
