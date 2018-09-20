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
<<<<<<< HEAD
    bool operator()(args &&...);
=======
    bool operator()(args &&...) const;
>>>>>>> 41b42cc5d45c544a427ed079af957ad4366ea9e6
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
<<<<<<< HEAD
    void operator()(args &&... );
=======
    void operator()(args &&... ) const;
>>>>>>> 41b42cc5d45c544a427ed079af957ad4366ea9e6
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
<<<<<<< HEAD
    void operator()(t1 &&, t2 &&);

    template <typename t>
    bool operator()(t &&, t &&);
=======
    void operator()(t1 &&, t2 &&) const;

    template <typename t>
    bool operator()(t &&, t &&) const;
>>>>>>> 41b42cc5d45c544a427ed079af957ad4366ea9e6
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
