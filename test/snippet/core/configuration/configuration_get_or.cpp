// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/core/configuration/pipeable_config_element.hpp>
#include <seqan3/core/debug_stream.hpp>

// Initial setup used in the actual example:
enum struct my_id : int
{
    bar_id,
    foo_id
};

struct bar : private seqan3::pipeable_config_element
{
public:
    float value{};

    bar() = default;
    bar(bar const &) = default;
    bar(bar &&) = default;
    bar & operator=(bar const &) = default;
    bar & operator=(bar &&) = default;
    ~bar() = default;

    bar(float v) : value{v}
    {}

    static constexpr my_id id{my_id::bar_id};
};

template <typename t>
class foo : private seqan3::pipeable_config_element
{
public:
    t value{};

    foo() = default;
    foo(foo const &) = default;
    foo(foo &&) = default;
    foo & operator=(foo const &) = default;
    foo & operator=(foo &&) = default;
    ~foo() = default;

    foo(t v) : value{std::move(v)}
    {}

    static constexpr my_id id{my_id::foo_id};
};

template <typename t>
foo(t) -> foo<t>;

int main()
{
    seqan3::configuration my_cfg{foo{1}};                                           // Only foo<int> is present.
    seqan3::debug_stream << my_cfg.get_or(foo{std::string{"hello"}}).value << '\n'; // finds foo<int> -> prints: 1
    seqan3::debug_stream << my_cfg.get_or(bar{2.4}).value << '\n';                  // bar not present -> prints: 2.4
}
