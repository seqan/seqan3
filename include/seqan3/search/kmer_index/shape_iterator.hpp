// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \author Vinzenz May <vinzenz.may AT fu-berlin.de>
 * \brief Provides the seqan3::shape_iterator.
 */

#pragma once

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/search/kmer_index/shape.hpp>
#include <seqan3/std/ranges>

namespace seqan3
{

/*!\brief Iterator for calculating hash values via a given seqan3::shape.
 * \tparam it_t Type of iterator on text. Must model std::ForwardIterator. Value type must model seqan3::Semialphabet.
*/
template <std::ForwardIterator it_t>
    //!\cond
    requires Semialphabet<value_type_t<it_t>>
    //!\endcond
class shape_iterator
{
private:
    //!\brief The alphabet type of the passed iterator.
    using alphabet_t = value_type_t<it_t>;

    //!\brief The hash value.
    size_t hash_value{0};

    //!\brief The seqan3::shape to use for hashing.
    shape const s;

    //!\brief The factor for the left most position of the hash value.
    size_t roll_factor{0};

    //!\brief True if `s` contains only 1 (ungapped), false otherwise.
    bool ungapped_shape{false};

    //!\brief Iterator pointing to begin of the underlying text.
    it_t const text_start;

    //!\brief Calculates a hash value by explicitly looking at each position.
    void hash_full()
    {
        text_right = text_left;
        hash_value = s[0] * to_rank(*(text_right));

        for (size_t i{1}; i < std::ranges::size(s); ++i)
        {
            std::advance(text_right, 1);

            hash_value *= alphabet_size_v<alphabet_t>;
            hash_value += s[i] * to_rank(*(text_right));
        }
    }

    //!\brief Calculates a hash value by using rolling hash.
    void hash_roll()
    {
        hash_value -= to_rank(*(text_left)) * roll_factor;

        std::advance(text_left,  1);
        std::advance(text_right, 1);

        hash_value *= alphabet_size_v<alphabet_t>;
        hash_value += to_rank(*(text_right));
    }

public:
    //!\brief Iterator to the leftmost position of the k-mer.
    it_t text_left;
    //!\brief Iterator to the rightmost position of the k-mer.
    it_t text_right;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    shape_iterator()                                   = default; //!< Defaulted.
    shape_iterator(shape_iterator const &)             = default; //!< Defaulted.
    shape_iterator(shape_iterator &&)                  = default; //!< Defaulted.
    shape_iterator & operator=(shape_iterator const &) = default; //!< Defaulted.
    shape_iterator & operator=(shape_iterator &&)      = default; //!< Defaulted.
    ~shape_iterator()                                  = default; //!< Defaulted.

    /*!\brief Construct from a given iterator on the text and a seqan3::shape.
    * /param[in] it_start Iterator pointing to the first position of the text.
    * /param[in] s_ seqan3::shape
    *
    * ### Complexity
    *
    * Linear in size of shape.
    */
    shape_iterator(it_t it_start, shape s_) : s(s_),
                                    text_start(it_start),
                                    text_left(it_start),
                                    text_right(it_start)
    {
        assert(std::ranges::size(s) > 0);
        ungapped_shape = std::all_of(std::ranges::begin(s),
                                     std::next(std::ranges::begin(s), std::ranges::size(s)),
                                     [](bool b) { return b; } );

        roll_factor = std::pow(alphabet_size_v<alphabet_t>, std::ranges::size(s) - 1);

        hash_full();
    };
    //!\}

    //!\name Comparison operators
    //!\{
    //\!brief Compare to iterator on text.
    inline bool operator==(it_t const & rhs) const noexcept
    {
        return text_right == rhs;
    }

    //\!brief Compare to another shape_iterator.
    inline bool operator==(shape_iterator const & rhs) const noexcept
    {
        return text_right == rhs.text_left;
    }

    //\!brief Compare to iterator on text.
    inline bool operator!=(it_t const & rhs) const noexcept
    {
        return !(text_right == rhs);
    }

    //\!brief Compare to another shape_iterator.
    inline bool operator!=(shape_iterator const & rhs) const noexcept
    {
        return !(text_right == rhs.text_left);
    }
    //!\}

    //!\brief Pre increment.
    shape_iterator & operator++()
    {
        if (ungapped_shape)
        {
            hash_roll();
        }
        else
        {
            std::advance(text_left,  1);
            hash_full();
        }
        return *this;
    }

    //!\brief Post increment.
    shape_iterator & operator++(int)
    {
        shape_iterator tmp{*this};
        if (ungapped_shape)
        {
            hash_roll();
        }
        else
        {
            std::advance(text_left,  1);
            hash_full();
        }
        return tmp;
    }

    //!\brief Calculate a hash value at a given position.
    shape_iterator & operator[](int offset)
        //!\cond
        requires std::RandomAccessIterator<it_t>
        //!\endcond
    {
        text_left = text_start + offset;

        hash_full();

        return *this;
    }


    //!\brief Return the hash value.
    size_t operator*()
    {
        return hash_value;
    }
};

}// namespace seqan3
