// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
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
// Author: Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
// ============================================================================

#pragma once

#include <algorithm>
#include <iostream>
#include <vector>
#include <array>
#include <list>
#include <forward_list>
#include <deque>
#include <string>

#include <seqan3/alphabet/concept.hpp>

namespace seqan3
{

/*! \file range.hpp
 * \ingroup alphabet
 * Free function overloads for containers over alphabets
 */

// ------------------------------------------------------------------
// ostream operator
// ------------------------------------------------------------------

//!\publicsection
//!@name outstream operator specializations for containers of our alphabets
//!@{
template <typename alphabet_type, size_t n>
    requires alphabet_concept<alphabet_type>
std::ostream& operator<<(std::ostream & os, std::array<alphabet_type, n> const & str)
{
    for (auto c : str)
        os << c;
    return os;
}

template <typename alphabet_type>
    requires alphabet_concept<alphabet_type>
std::ostream& operator<<(std::ostream & os, std::forward_list<alphabet_type> const & str)
{
    for (auto c : str)
        os << c;
    return os;
}

template <typename alphabet_type>
    requires alphabet_concept<alphabet_type>
std::ostream& operator<<(std::ostream & os, std::list<alphabet_type> const & str)
{
    for (auto c : str)
        os << c;
    return os;
}

template <typename alphabet_type>
    requires alphabet_concept<alphabet_type>
std::ostream& operator<<(std::ostream & os, std::deque<alphabet_type> const & str)
{
    for (auto c : str)
        os << c;
    return os;
}

template <typename alphabet_type>
    requires alphabet_concept<alphabet_type>
std::ostream& operator<<(std::ostream & os, std::vector<alphabet_type> const & str)
{
    for (auto c : str)
        os << c;
    return os;
}

template <typename alphabet_type>
    requires alphabet_concept<alphabet_type>
std::ostream& operator<<(std::ostream & os, std::basic_string<alphabet_type> const & str)
{
    for (auto c : str)
        os << c;
    return os;
}
//!@}
//TODO serialization

}
