// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include "gtest/gtest.h"

#include <array>

#include <seqan3/core/algorithm/configuration_utility.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>

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
                            static_cast<uint8_t>(::test_algo_id::SIZE)> compatibility_table<::test_algo_id>
{{
    {0, 1, 1, 1},
    {1, 0, 1, 0},
    {1, 1, 0, 1},
    {1, 0, 1, 0}
}};

}  // namespace seqan3::detail

struct bar : public seqan3::pipeable_config_element<bar, int>
{

    /*!\name Constructors, destructor and assignment
     * \{
     */
    bar() = default; //!< Defaulted.
    bar(bar const &) = default; //!< Defaulted.
    bar(bar &&) = default; //!< Defaulted.
    bar & operator=(bar const &) = default; //!< Defaulted.
    bar & operator=(bar &&) = default; //!< Defaulted.
    ~bar() = default; //!< Defaulted.

    //!\brief Construct from base type.
    constexpr bar(int const & i) : seqan3::pipeable_config_element<bar, int>(i) {}
    //!\}

    static constexpr test_algo_id id{test_algo_id::bar_id};
};

struct bax : public seqan3::pipeable_config_element<bax, float>
{
    /*!\name Constructors, destructor and assignment
     * \{
     */
    bax() = default; //!< Defaulted.
    bax(bax const &) = default; //!< Defaulted.
    bax(bax &&) = default; //!< Defaulted.
    bax & operator=(bax const &) = default; //!< Defaulted.
    bax & operator=(bax &&) = default; //!< Defaulted.
    ~bax() = default; //!< Defaulted.

    //!\brief Construct from base type.
    constexpr bax(float const & f) : seqan3::pipeable_config_element<bax, float>(f) {}
    //!\}

    static constexpr test_algo_id id{test_algo_id::bax_id};
};

struct foo : public seqan3::pipeable_config_element<foo, std::string>
{
    /*!\name Constructors, destructor and assignment
     * \{
     */
    foo() = default; //!< Defaulted.
    foo(foo const &) = default; //!< Defaulted.
    foo(foo &&) = default; //!< Defaulted.
    foo & operator=(foo const &) = default; //!< Defaulted.
    foo & operator=(foo &&) = default; //!< Defaulted.
    ~foo() = default; //!< Defaulted.

    //!\brief Construct from base type.
    constexpr foo(std::string const & s) : seqan3::pipeable_config_element<foo, std::string>(s) {}
    //!\}

    static constexpr test_algo_id id{test_algo_id::foo_id};
};

template <typename t = std::vector<int>>
struct foobar : public seqan3::pipeable_config_element<foobar<t>, t>
{
    /*!\name Constructors, destructor and assignment
     * \{
     */
    foobar() = default; //!< Defaulted.
    foobar(foobar const &) = default; //!< Defaulted.
    foobar(foobar &&) = default; //!< Defaulted.
    foobar & operator=(foobar const &) = default; //!< Defaulted.
    foobar & operator=(foobar &&) = default; //!< Defaulted.
    ~foobar() = default; //!< Defaulted.

    //!\brief Construct from base type.
    constexpr foobar(t const & e) : seqan3::pipeable_config_element<foobar<t>, t>(e) {}
    //!\}

    static constexpr test_algo_id id{test_algo_id::foobar_id};
};
