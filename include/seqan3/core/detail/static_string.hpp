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

/*!\file
 * \ingroup core
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \brief A static string implementation to manipulate string literals at compile time.
 */

#include <type_traits>
#include <array>

namespace seqan3
{

template <std::size_t N>
class static_string
{
public:

    static_string() = delete;
    static_string(static_string const &) = default;
    static_string(static_string &&) = default;
    static_string & operator=(static_string const &) = default;
    static_string & operator=(static_string &&) = default;

    // Construction from two static_strings.
    template <std::size_t N1>
    constexpr static_string(static_string<N1> const & lhs,
                            static_string<N - N1> const & rhs) noexcept
    {
        for (unsigned i = 0; i < N1; ++i)
            _lit[i] = lhs[i];
        for (unsigned i = N1; i < N + 1; ++i)
            _lit[i] = rhs[i - N1];
    }

    // Construction from literal.
    constexpr static_string(const char (&lit)[N + 1]) noexcept
    {
        // static_assert(lit[N] == '\0'); TODO(rrahn): Fix me
        for (unsigned i = 0; i < N + 1; ++i)
            _lit[i] = lit[i];
    }

    // Construction from char.
    constexpr static_string(const char c) noexcept
    {
        _lit[0] = c;
        _lit[1] = '\0';
    }

    // access via subscript
    constexpr char operator[](std::size_t n) const noexcept
    {
        return _lit[n];
    }

    // Concatenation
    template <std::size_t N2>
    constexpr static_string<N + N2> operator+(static_string<N2> const & rhs) noexcept
    {
        return static_string<N + N2>{*this, rhs};
    }

    constexpr std::size_t size() noexcept
    {
        return N;
    }

    std::string string()
    {
        return std::string{ _lit.cbegin(), _lit.cend() - 1};
    }

    constexpr const char * c_str() noexcept
    {
        return _lit.data();
    }

protected:
  std::array<char, N + 1> _lit{};
};

// User defined deduction guide for string literals
template <std::size_t N>
static_string(const char (&)[N]) -> static_string<N - 1>;

// User defined deduction guide for chars
static_string(const char) -> static_string<1>;

}  // namespace seqan3
