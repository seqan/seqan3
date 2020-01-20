// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

//! \brief Helper struct for testing core concepts.
struct type_a
{};

//! \brief Helper struct for testing core concepts.
struct type_b : type_a
{
    type_b(type_b const &) = delete;
    type_b(type_b &&) noexcept = default;

    type_b & operator=(type_b const &) = delete;
    type_b & operator=(type_b &&) noexcept = default;

    template <typename ...args>
    bool operator()(args &&...) const;
};

//! \brief Helper struct for testing core concepts.
struct type_c
{
    type_c() = default;
    type_c(type_b const &)
    {}
    explicit type_c(type_a const &)
    {}

    template <typename ...args>
    void operator()(args &&... ) const;
};

//! \brief Helper struct for testing core concepts.
struct type_d: type_b
{
    type_d() = delete;

    type_d(type_d const &) = delete;
    type_d(type_d &&) = delete;

    type_d & operator=(type_d const &) = delete;
    type_d & operator=(type_d &&) = delete;
    ~type_d() = delete;

    template <typename t1, typename t2>
    void operator()(t1 &&, t2 &&) const;

    template <typename t>
    bool operator()(t &&, t &&) const;
};

// Operator overloads for testing core concepts.
bool operator==(type_a const & , type_b const &);
bool operator==(type_b const & , type_a const &);
bool operator==(type_b const & , type_b const &);
bool operator==(type_d const & , type_b const &);
bool operator==(type_b const & , type_d const &);
bool operator==(type_d const & , type_d const &);
bool operator==(type_c const & , type_c const &);

bool operator!=(type_a const & , type_b const &);
bool operator!=(type_b const & , type_a const &);
bool operator!=(type_b const & , type_b const &);
bool operator!=(type_d const & , type_b const &);
bool operator!=(type_b const & , type_d const &);
bool operator!=(type_d const & , type_d const &);
bool operator!=(type_c const & , type_c const &);

bool operator<(type_a const & , type_a const &);
bool operator<(type_a const & , type_b const &);
bool operator<(type_b const & , type_b const &);
bool operator<(type_b const & , type_a const &);
bool operator>(type_a const & , type_a const &);
bool operator>(type_a const & , type_b const &);
bool operator>(type_b const & , type_b const &);
bool operator>(type_b const & , type_a const &);
bool operator<=(type_a const & , type_a const &);
bool operator<=(type_a const & , type_b const &);
bool operator<=(type_b const & , type_b const &);
bool operator<=(type_b const & , type_a const &);
bool operator>=(type_a const & , type_a const &);
bool operator>=(type_a const & , type_b const &);
bool operator>=(type_b const & , type_b const &);
bool operator>=(type_b const & , type_a const &);

bool operator<(type_d const & , type_d const &);
bool operator<(type_d const & , type_b const &);
bool operator<(type_b const & , type_d const &);
bool operator>(type_d const & , type_d const &);
bool operator>(type_d const & , type_b const &);
bool operator>(type_b const & , type_d const &);
bool operator<=(type_d const & , type_d const &);
bool operator<=(type_d const & , type_b const &);
bool operator<=(type_b const & , type_d const &);
bool operator>=(type_d const & , type_d const &);
bool operator>=(type_d const & , type_b const &);
bool operator>=(type_b const & , type_d const &);
