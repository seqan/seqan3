// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <seqan3/std/memory>

#include <gtest/gtest.h>

#include <seqan3/test/expect_same_type.hpp>

TEST(to_address, pointer)
{
    int value{};
    int * value_ptr = &value;

    EXPECT_EQ(value_ptr, std::to_address(value_ptr));
}

TEST(to_address, pointer_traits)
{
    std::shared_ptr<int> value_ptr = std::make_shared<int>(5);

    EXPECT_SAME_TYPE(std::pointer_traits<std::shared_ptr<int>>::pointer,
                     std::shared_ptr<int>);

    EXPECT_SAME_TYPE(std::pointer_traits<std::shared_ptr<int>>::element_type,
                     int);

    EXPECT_EQ(*value_ptr, 5);
    EXPECT_EQ(value_ptr.get(), std::to_address(value_ptr));
}

struct fancy_ptr_t
{
    using element_type = int; // needed to enable std::pointer_traits

    int * operator->()
    {
        return &value;
    }

    int const * operator->() const
    {
        return &value;
    }

    int value{5};
};

TEST(to_address, pointer_traits_and_member_arrow_operator)
{
    fancy_ptr_t fancy_ptr{};
    int * value_ptr = &fancy_ptr.value;

    EXPECT_SAME_TYPE(std::pointer_traits<fancy_ptr_t>::pointer,
                     fancy_ptr_t);

    EXPECT_SAME_TYPE(std::pointer_traits<fancy_ptr_t>::element_type,
                     int);

    EXPECT_EQ(value_ptr, std::to_address(fancy_ptr));
}

struct fancy_ptr2_t
{
    int * value_ptr;
};

namespace std
{
template <>
struct pointer_traits<fancy_ptr2_t>
{
    using pointer = fancy_ptr2_t;
    using element_type = int;
    using difference_type = std::ptrdiff_t;

    static element_type * to_address(pointer p) noexcept
    {
        return p.value_ptr;
    }
};
} // namespace std

TEST(to_address, pointer_traits_to_address)
{
    int value{5};
    fancy_ptr2_t fancy_ptr{&value};

    EXPECT_SAME_TYPE(std::pointer_traits<fancy_ptr2_t>::pointer,
                     fancy_ptr2_t);

    EXPECT_SAME_TYPE(std::pointer_traits<fancy_ptr2_t>::element_type,
                     int);

    EXPECT_EQ(fancy_ptr.value_ptr, std::pointer_traits<fancy_ptr2_t>::to_address(fancy_ptr));
    EXPECT_EQ(fancy_ptr.value_ptr, std::to_address(fancy_ptr));
}
