// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/core/configuration/pipeable_config_element.hpp>

enum struct my_id : int
{
    bar_id,
    foo_id
};

class bar : private seqan3::pipeable_config_element
{
public:
    bar() = default;
    bar(bar const &) = default;
    bar(bar &&) = default;
    bar & operator=(bar const &) = default;
    bar & operator=(bar &&) = default;
    ~bar() = default;

    static constexpr my_id id{my_id::bar_id};
};

template <typename t>
struct foo : private seqan3::pipeable_config_element
{
public:
    foo() = default;
    foo(foo const &) = default;
    foo(foo &&) = default;
    foo & operator=(foo const &) = default;
    foo & operator=(foo &&) = default;
    ~foo() = default;

    static constexpr my_id id{my_id::foo_id};
};

template <typename t>
foo(t) -> foo<t>;
