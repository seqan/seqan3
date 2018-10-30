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

#include <sstream>
#include <vector>

#include <gtest/gtest.h>

#include <seqan3/range/shortcuts.hpp>
#include <seqan3/range/detail/inherited_iterator_base.hpp>

/* This class is extensively tested by the many views that use it, e.g. view::take_line */

using namespace seqan3;

//! [inherited_iterator_base def]

class skip_odd_numbers_it : public detail::inherited_iterator_base<skip_odd_numbers_it, std::vector<int>::iterator>
{
private:
    using base_base_t = std::vector<int>::iterator;
    using base_t      = detail::inherited_iterator_base<skip_odd_numbers_it, std::vector<int>::iterator>;

public:
    skip_odd_numbers_it() = default;
    constexpr skip_odd_numbers_it(skip_odd_numbers_it const & rhs) = default;
    constexpr skip_odd_numbers_it(skip_odd_numbers_it && rhs) = default;
    constexpr skip_odd_numbers_it & operator=(skip_odd_numbers_it const & rhs) = default;
    constexpr skip_odd_numbers_it & operator=(skip_odd_numbers_it && rhs) = default;
    ~skip_odd_numbers_it() = default;

    skip_odd_numbers_it(base_base_t it) : base_t{it} {}

    using difference_type       = typename std::iterator_traits<base_base_t>::difference_type;
    using value_type            = typename std::iterator_traits<base_base_t>::value_type;
    using reference             = typename std::iterator_traits<base_base_t>::reference;
    using pointer               = typename std::iterator_traits<base_base_t>::pointer;
    using iterator_category     = typename std::iterator_traits<base_base_t>::iterator_category;

    skip_odd_numbers_it & operator++()
    {
        *this = ++static_cast<base_base_t>(*this);

        if ((**this) % 2 != 0)
            *this = ++static_cast<base_base_t>(*this);

        return *this;
    }

    skip_odd_numbers_it operator++(int)
    {
        skip_odd_numbers_it cpy{*this};
        ++(*this);
        return cpy;
    }
};
//! [inherited_iterator_base def]

TEST(inherited_iterator_base, minimal)
{
    //! [inherited_iterator_base desired]
    std::vector<int> vec{0,1,2,3,4,5,6,7,8,9};

    skip_odd_numbers_it it = begin(vec);

    EXPECT_EQ(*it, 0);
    ++it;
    EXPECT_EQ(*it, 2);
    ++it;
    EXPECT_EQ(*it, 4);
    ++it;
    EXPECT_EQ(*it, 6);
    ++it;
    EXPECT_EQ(*it, 8);
    //! [inherited_iterator_base desired]
}

TEST(inherited_iterator_base, concept_check)
{
    EXPECT_TRUE(std::RandomAccessIterator<skip_odd_numbers_it>);
}
