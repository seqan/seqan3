// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <sstream>
#include <vector>

#include <seqan3/core/range/detail/inherited_iterator_base.hpp>

/* This class is extensively tested by the many views that use it, e.g. views::take_line */

//! [inherited_iterator_base def]

class skip_odd_numbers_it :
    public seqan3::detail::inherited_iterator_base<skip_odd_numbers_it, std::vector<int>::iterator>
{
private:
    using base_base_t = std::vector<int>::iterator;
    using base_t = seqan3::detail::inherited_iterator_base<skip_odd_numbers_it, std::vector<int>::iterator>;

public:
    skip_odd_numbers_it() = default;
    skip_odd_numbers_it(skip_odd_numbers_it const & rhs) = default;
    skip_odd_numbers_it(skip_odd_numbers_it && rhs) = default;
    skip_odd_numbers_it & operator=(skip_odd_numbers_it const & rhs) = default;
    skip_odd_numbers_it & operator=(skip_odd_numbers_it && rhs) = default;
    ~skip_odd_numbers_it() = default;

    skip_odd_numbers_it(base_base_t it) : base_t{it}
    {}

    using difference_type = typename std::iterator_traits<base_base_t>::difference_type;
    using value_type = typename std::iterator_traits<base_base_t>::value_type;
    using reference = typename std::iterator_traits<base_base_t>::reference;
    using pointer = typename std::iterator_traits<base_base_t>::pointer;
    using iterator_category = typename std::iterator_traits<base_base_t>::iterator_category;

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
    std::vector<int> vec{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

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
    EXPECT_TRUE(std::random_access_iterator<skip_odd_numbers_it>);
}
