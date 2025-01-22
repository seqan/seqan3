// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <array>

#include <seqan3/core/configuration/detail/concept.hpp>
#include <seqan3/core/configuration/pipeable_config_element.hpp>

enum class test_algo_id : uint8_t
{
    bar_id = 0,
    bax_id = 1,
    foo_id = 2,
    foobar_id = 3,
    SIZE = 4
};

namespace seqan3::detail
{

template <>
inline constexpr std::array<std::array<bool, static_cast<uint8_t>(::test_algo_id::SIZE)>,
                            static_cast<uint8_t>(::test_algo_id::SIZE)>
    compatibility_table<::test_algo_id>{{{0, 1, 1, 1}, {1, 0, 1, 0}, {1, 1, 0, 1}, {1, 0, 1, 0}}};

} // namespace seqan3::detail

class bar : private seqan3::pipeable_config_element
{
public:
    int value{};

    constexpr bar() = default;
    constexpr bar(bar const &) = default;
    constexpr bar(bar &&) = default;
    constexpr bar & operator=(bar const &) = default;
    constexpr bar & operator=(bar &&) = default;
    ~bar() = default;

    constexpr bar(int v) : value{v}
    {}

    static constexpr test_algo_id id{test_algo_id::bar_id};
};

class bax : private seqan3::pipeable_config_element
{
public:
    float value{};

    constexpr bax() = default;
    constexpr bax(bax const &) = default;
    constexpr bax(bax &&) = default;
    constexpr bax & operator=(bax const &) = default;
    constexpr bax & operator=(bax &&) = default;
    ~bax() = default;

    constexpr bax(float v) : value{v}
    {}

    static constexpr test_algo_id id{test_algo_id::bax_id};
};

class foo : private seqan3::pipeable_config_element
{
public:
    std::string value{};

    foo() = default;
    foo(foo const &) = default;
    foo(foo &&) = default;
    foo & operator=(foo const &) = default;
    foo & operator=(foo &&) = default;
    ~foo() = default;

    foo(std::string v) : value{v}
    {}

    static constexpr test_algo_id id{test_algo_id::foo_id};
};

template <typename t = std::vector<int>>
class foobar : private seqan3::pipeable_config_element
{
public:
    t value{};

    constexpr foobar() = default;
    constexpr foobar(foobar const &) = default;
    constexpr foobar(foobar &&) = default;
    constexpr foobar & operator=(foobar const &) = default;
    constexpr foobar & operator=(foobar &&) = default;
    ~foobar() = default;

    constexpr foobar(t v) : value{v}
    {}

    static constexpr test_algo_id id{test_algo_id::foobar_id};
};
