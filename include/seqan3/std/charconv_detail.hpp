// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

// ============================================================================
// LLVM Release License
// ============================================================================
// University of Illinois/NCSA
// Open Source License
//
// Copyright (c) 2003-2018 University of Illinois at Urbana-Champaign.
// All rights reserved.
//
// Developed by:
//     LLVM Team
//     University of Illinois at Urbana-Champaign
//     http://llvm.org
//
// Permission is hereby granted, free of charge, to any person obtaining a copy of
// this software and associated documentation files (the "Software"), to deal with
// the Software without restriction, including without limitation the rights to
// use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
// of the Software, and to permit persons to whom the Software is furnished to do
// so, subject to the following conditions:
//
//     * Redistributions of source code must retain the above copyright notice,
//       this list of conditions and the following disclaimers.
//     * Redistributions in binary form must reproduce the above copyright notice,
//       this list of conditions and the following disclaimers in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the names of the LLVM Team, University of Illinois at
//       Urbana-Champaign, nor the names of its contributors may be used to
//       endorse or promote products derived from this Software without specific
//       prior written permission.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
// FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
// CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS WITH THE
// SOFTWARE.
// ============================================================================

/*!\file
 * \brief Provides LLVM's implementation of std::from_chars and std::to_chars
 * for integral types and a custom solution for floating points types.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <type_traits>

#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/io/stream/parse_condition.hpp>
#include <seqan3/std/concepts>

namespace seqan3::detail
{

//!\cond
// implementation detail taken from LLVM
static constexpr uint64_t pow10_64[] =
{
    UINT64_C(0),
    UINT64_C(10),
    UINT64_C(100),
    UINT64_C(1000),
    UINT64_C(10000),
    UINT64_C(100000),
    UINT64_C(1000000),
    UINT64_C(10000000),
    UINT64_C(100000000),
    UINT64_C(1000000000),
    UINT64_C(10000000000),
    UINT64_C(100000000000),
    UINT64_C(1000000000000),
    UINT64_C(10000000000000),
    UINT64_C(100000000000000),
    UINT64_C(1000000000000000),
    UINT64_C(10000000000000000),
    UINT64_C(100000000000000000),
    UINT64_C(1000000000000000000),
    UINT64_C(10000000000000000000),
};

static constexpr uint32_t pow10_32[] =
{
    UINT32_C(0),          UINT32_C(10),       UINT32_C(100),
    UINT32_C(1000),       UINT32_C(10000),    UINT32_C(100000),
    UINT32_C(1000000),    UINT32_C(10000000), UINT32_C(100000000),
    UINT32_C(1000000000),
};

static constexpr char cDigitsLut[200] = {
    '0', '0', '0', '1', '0', '2', '0', '3', '0', '4', '0', '5', '0', '6', '0',
    '7', '0', '8', '0', '9', '1', '0', '1', '1', '1', '2', '1', '3', '1', '4',
    '1', '5', '1', '6', '1', '7', '1', '8', '1', '9', '2', '0', '2', '1', '2',
    '2', '2', '3', '2', '4', '2', '5', '2', '6', '2', '7', '2', '8', '2', '9',
    '3', '0', '3', '1', '3', '2', '3', '3', '3', '4', '3', '5', '3', '6', '3',
    '7', '3', '8', '3', '9', '4', '0', '4', '1', '4', '2', '4', '3', '4', '4',
    '4', '5', '4', '6', '4', '7', '4', '8', '4', '9', '5', '0', '5', '1', '5',
    '2', '5', '3', '5', '4', '5', '5', '5', '6', '5', '7', '5', '8', '5', '9',
    '6', '0', '6', '1', '6', '2', '6', '3', '6', '4', '6', '5', '6', '6', '6',
    '7', '6', '8', '6', '9', '7', '0', '7', '1', '7', '2', '7', '3', '7', '4',
    '7', '5', '7', '6', '7', '7', '7', '8', '7', '9', '8', '0', '8', '1', '8',
    '2', '8', '3', '8', '4', '8', '5', '8', '6', '8', '7', '8', '8', '8', '9',
    '9', '0', '9', '1', '9', '2', '9', '3', '9', '4', '9', '5', '9', '6', '9',
    '7', '9', '8', '9', '9'};

template <typename T>
inline char* append1(char* buffer, T i) noexcept
{
    *buffer = '0' + static_cast<char>(i);
    return buffer + 1;
}

template <typename T>
inline char* append2(char* buffer, T i) noexcept
{
    memcpy(buffer, &cDigitsLut[(i)*2], 2);
    return buffer + 2;
}

template <typename T>
inline char* append3(char* buffer, T i) noexcept
{
    return append2(append1(buffer, (i) / 100), (i) % 100);
}

template <typename T>
inline char* append4(char* buffer, T i) noexcept
{
    return append2(append2(buffer, (i) / 100), (i) % 100);
}

inline char* u32toa(uint32_t value, char* buffer) noexcept
{
    if (value < 10000)
    {
        if (value < 100)
        {
            if (value < 10)
                buffer = append1(buffer, value);
            else
                buffer = append2(buffer, value);
        }
        else
        {
            if (value < 1000)
                buffer = append3(buffer, value);
            else
                buffer = append4(buffer, value);
        }
    }
    else if (value < 100000000)
    {
        // value = bbbbcccc
        const uint32_t b = value / 10000;
        const uint32_t c = value % 10000;

        if (value < 1000000)
        {
            if (value < 100000)
                buffer = append1(buffer, b);
            else
                buffer = append2(buffer, b);
        }
        else
        {
            if (value < 10000000)
                buffer = append3(buffer, b);
            else
                buffer = append4(buffer, b);
        }

        buffer = append4(buffer, c);
    }
    else
    {
        // value = aabbbbcccc in decimal
        const uint32_t a = value / 100000000;  // 1 to 42
        value %= 100000000;

        if (a < 10)
            buffer = append1(buffer, a);
        else
            buffer = append2(buffer, a);

        buffer = append4(buffer, value / 10000);
        buffer = append4(buffer, value % 10000);
    }

    return buffer;
}

inline char* u64toa(uint64_t value, char* buffer) noexcept
{
    if (value < 100000000)
    {
        uint32_t v = static_cast<uint32_t>(value);
        if (v < 10000)
        {
            if (v < 100)
            {
                if (v < 10)
                    buffer = append1(buffer, v);
                else
                    buffer = append2(buffer, v);
            }
            else
            {
                if (v < 1000)
                    buffer = append3(buffer, v);
                else
                    buffer = append4(buffer, v);
            }
        }
        else
        {
            // value = bbbbcccc
            const uint32_t b = v / 10000;
            const uint32_t c = v % 10000;

            if (v < 1000000)
            {
                if (v < 100000)
                    buffer = append1(buffer, b);
                else
                    buffer = append2(buffer, b);
            }
            else
            {
                if (v < 10000000)
                    buffer = append3(buffer, b);
                else
                    buffer = append4(buffer, b);
            }

            buffer = append4(buffer, c);
        }
    }
    else if (value < 10000000000000000)
    {
        const uint32_t v0 = static_cast<uint32_t>(value / 100000000);
        const uint32_t v1 = static_cast<uint32_t>(value % 100000000);

        const uint32_t b0 = v0 / 10000;
        const uint32_t c0 = v0 % 10000;

        if (v0 < 1000000)
        {
            if (v0 < 100000)
                buffer = append1(buffer, b0);
            else
                buffer = append2(buffer, b0);
        }
        else
        {
            if (v0 < 10000000)
                buffer = append3(buffer, b0);
            else
                buffer = append4(buffer, b0);
        }

        buffer = append4(buffer, c0);
        buffer = append4(buffer, v1 / 10000);
        buffer = append4(buffer, v1 % 10000);
    }
    else
    {
        const uint32_t a =
            static_cast<uint32_t>(value / 10000000000000000);  // 1 to 1844
        value %= 10000000000000000;

        if (a < 100)
        {
            if (a < 10)
                buffer = append1(buffer, a);
            else
                buffer = append2(buffer, a);
        }
        else
        {
            if (a < 1000)
                buffer = append3(buffer, a);
            else
                buffer = append4(buffer, a);
        }

        const uint32_t v0 = static_cast<uint32_t>(value / 100000000);
        const uint32_t v1 = static_cast<uint32_t>(value % 100000000);
        buffer = append4(buffer, v0 / 10000);
        buffer = append4(buffer, v0 % 10000);
        buffer = append4(buffer, v1 / 10000);
        buffer = append4(buffer, v1 % 10000);
    }

    return buffer;
}

template <typename value_type, typename = void>
struct traits_base
{
    using type = uint64_t;

    static int width(value_type v) noexcept
    {
        auto t = (64 - __builtin_clzll(v | 1)) * 1233 >> 12;
        return t - (v < pow10_64[t]) + 1;
    }

    static char* convert(value_type v, char* p) noexcept
    {
        return u64toa(v, p);
    }

    static auto& pow() { return pow10_64; }
};

template <typename value_type>
struct traits_base<value_type, decltype(void(uint32_t{std::declval<value_type>()}))>
{
    using type = uint32_t;

    static int width(value_type v) noexcept
    {
        auto t = (32 - __builtin_clz(v | 1)) * 1233 >> 12;
        return t - (v < pow10_32[t]) + 1;
    }

    static char* convert(value_type v, char* p) noexcept
    {
        return u32toa(v, p);
    }

    static auto& pow() noexcept { return pow10_32; }
};

template <typename value_type>
inline bool mul_overflowed(unsigned char a, value_type b, unsigned char& r) noexcept
{
    auto c = a * b;
    r = c;
    return c > (std::numeric_limits<unsigned char>::max)();
}

template <typename value_type>
inline bool mul_overflowed(unsigned short a, value_type b, unsigned short& r) noexcept
{
    auto c = a * b;
    r = c;
    return c > (std::numeric_limits<unsigned short>::max)();
}

template <typename value_type>
inline bool mul_overflowed(value_type a, value_type b, value_type & r) noexcept
{
    static_assert(std::is_unsigned<value_type>::value, "");
    return __builtin_mul_overflow(a, b, &r);
}

template <typename value_type, typename _Up>
inline bool mul_overflowed(value_type a, _Up b, value_type & r) noexcept
{
    return mul_overflowed(a, static_cast<value_type>(b), r);
}

template <typename value_type>
struct traits : traits_base<value_type>
{
    static constexpr int digits = std::numeric_limits<value_type>::digits10 + 1;
    using traits_base<value_type>::pow;
    using typename traits_base<value_type>::type;

    // precondition: at least one non-zero character available
    static char const*
    read(char const* p, char const* ep, type& a, type& b) noexcept
    {
        type cprod[digits];
        int j = digits - 1;
        int i = digits;
        do
        {
            if (!('0' <= *p && *p <= '9'))
                break;
            cprod[--i] = *p++ - '0';
        } while (p != ep && i != 0);

        a = inner_product(cprod + i + 1, cprod + j, pow() + 1,
                              cprod[i]);
        if (mul_overflowed(cprod[j], pow()[j - i], b))
            --p;
        return p;
    }

    template <typename _It1, typename _It2, class _Up>
    static _Up
    inner_product(_It1 first1, _It1 last1, _It2 first2, _Up init) noexcept
    {
        for (; first1 < last1; ++first1, ++first2)
            init = init + *first1 * *first2;
        return init;
    }
};

template <typename value_type>
inline value_type complement(value_type x) noexcept
{
    static_assert(std::UnsignedIntegral<value_type>, "cast to unsigned first");
    return value_type(~x + 1);
}

template <typename value_type>
inline auto to_unsigned(value_type x) noexcept
{
    return static_cast<std::make_unsigned_t<value_type>>(x);
}

template <typename value_type>
inline std::to_chars_result to_chars_itoa(char* first, char* last, value_type value, std::false_type) noexcept
{
    using tx = traits<value_type>;
    auto diff = last - first;

    if (tx::digits <= diff || tx::width(value) <= diff)
        return {tx::convert(value, first), {}};
    else
        return {last, std::errc::value_too_large};
}

template <typename value_type>
inline std::to_chars_result to_chars_itoa(char* first, char* last, value_type value, std::true_type) noexcept
{
    auto x = to_unsigned(value);
    if (value < 0 && first != last)
    {
        *first++ = '-';
        x = complement(x);
    }

    return to_chars_itoa(first, last, x, std::false_type());
}

template <typename value_type>
inline std::to_chars_result to_chars_integral(char* first, char* last, value_type value, int base,
                    std::true_type) noexcept
{
    auto x = to_unsigned(value);
    if (value < 0 && first != last)
    {
        *first++ = '-';
        x = complement(x);
    }

    return to_chars_integral(first, last, x, base, std::false_type());
}

template <typename value_type>
inline std::to_chars_result to_chars_integral(char* first, char* last, value_type value, int base,
                    std::false_type) noexcept
{
    if (base == 10)
        return to_chars_itoa(first, last, value, std::false_type());

    auto p = last;
    while (p != first)
    {
        auto c = value % base;
        value /= base;
        *--p = "0123456789abcdefghijklmnopqrstuvwxyz"[c];
        if (value == 0)
            break;
    }

    auto len = last - p;
    if (value != 0 || !len)
        return {last, std::errc::value_too_large};
    else
    {
        memmove(first, p, len);
        return {first + len, {}};
    }
}

template <typename _It, typename value_type, typename _Fn, typename... _Ts>
inline std::from_chars_result sign_combinator(_It first, _It last, value_type & value, _Fn f, _Ts... args) noexcept
{
    using tl = std::numeric_limits<value_type>;
    decltype(to_unsigned(value)) x;

    bool neg = (first != last && *first == '-');
    auto r = f(neg ? first + 1 : first, last, x, args...);

    switch (r.ec)
    {
    case std::errc::invalid_argument:
        return {first, r.ec};
    case std::errc::result_out_of_range:
        return r;
    default:
        break;
    }

    if (neg)
    {
        if (x <= complement(to_unsigned(tl::min())))
        {
            x = complement(x);
            memcpy(&value, &x, sizeof(x));
            return r;
        }
    }
    else
    {
        if (x <= to_unsigned(tl::max()))
        {
            value = x;
            return r;
        }
    }

    return {r.ptr, std::errc::result_out_of_range};
}

template <typename value_type>
inline bool in_pattern(value_type c) noexcept
{
    return '0' <= c && c <= '9';
}

struct in_pattern_result
{
    bool ok;
    int val;

    explicit operator bool() const noexcept { return ok; }
};

template <typename value_type>
inline in_pattern_result in_pattern(value_type c, int base) noexcept
{
    if (base <= 10)
        return {'0' <= c && c < '0' + base, c - '0'};
    else if (in_pattern(c))
        return {true, c - '0'};
    else if ('a' <= c && c < 'a' + base - 10)
        return {true, c - 'a' + 10};
    else
        return {'A' <= c && c < 'A' + base - 10, c - 'A' + 10};
}

template <typename _It, typename value_type, typename _Fn, typename... _Ts>
inline std::from_chars_result subject_seq_combinator(_It first, _It last, value_type & value, _Fn f, _Ts... args) noexcept
{
    auto find_non_zero = [](_It first, _It last)
    {
        for (; first != last; ++first)
            if (*first != '0')
                break;
        return first;
    };

    auto p = find_non_zero(first, last);
    if (p == last || !in_pattern(*p, args...))
    {
        if (p == first)
            return {first, std::errc::invalid_argument};
        else
        {
            value = 0;
            return {p, {}};
        }
    }

    auto r = f(p, last, value, args...);
    if (r.ec == std::errc::result_out_of_range)
    {
        for (; r.ptr != last; ++r.ptr)
        {
            if (!in_pattern(*r.ptr, args...))
                break;
        }
    }

    return r;
}

template <typename value_type, std::enable_if_t<std::is_unsigned<value_type>::value, int> = 0>
inline std::from_chars_result
from_chars_atoi(char const * first, char const * last, value_type & value) noexcept
{
    using tx = traits<value_type>;
    using output_type = typename tx::type;

    return subject_seq_combinator(first, last, value,
                                  [](char const * first, char const * last, value_type & value) -> std::from_chars_result
                                  {
                                      output_type a, b;
                                      auto p = tx::read(first, last, a, b);
                                      if (p == last || !in_pattern(*p))
                                      {
                                          output_type m = (std::numeric_limits<value_type>::max)();
                                          if (m >= a && m - a >= b)
                                          {
                                              value = a + b;
                                              return {p, {}};
                                          }
                                      }
                                      return {p, std::errc::result_out_of_range};
                                  });
}

template <std::SignedIntegral value_type>
inline std::from_chars_result from_chars_atoi(char const * first, char const * last, value_type & value) noexcept
{
    using t = decltype(to_unsigned(value));
    return sign_combinator(first, last, value, from_chars_atoi<t>);
}

template <std::UnsignedIntegral value_type>
inline std::from_chars_result from_chars_integral(char const * first, char const * last, value_type & value, int base) noexcept
{
    if (base == 10)
        return from_chars_atoi(first, last, value);

    return subject_seq_combinator(first, last, value,
                                  [] (char const * p, char const * last, value_type & value, int base) -> std::from_chars_result
                                  {
                                      using tl = std::numeric_limits<value_type>;
                                      auto digits = tl::digits / log2f(float(base));
                                      value_type a = in_pattern(*p++, base).val, b = 0;

                                      for (int i = 1; p != last; ++i, ++p)
                                      {
                                          if (auto c = in_pattern(*p, base))
                                          {
                                              if (i < digits - 1)
                                                  a = a * base + c.val;
                                              else
                                              {
                                                  if (!mul_overflowed(a, base, a))
                                                      ++p;
                                                  b = c.val;
                                                  break;
                                              }
                                          }
                                          else
                                              break;
                                      }

                                      if (p == last || !in_pattern(*p, base))
                                      {
                                          if ((tl::max)() - a >= b)
                                          {
                                              value = a + b;
                                              return {p, {}};
                                          }
                                      }
                                      return {p, std::errc::result_out_of_range};
                                  },
                                  base);
}

template <std::SignedIntegral value_type>
inline std::from_chars_result from_chars_integral(char const * first, char const * last, value_type & value, int base) noexcept
{
    using t = decltype(to_unsigned(value));
    return sign_combinator(first, last, value, from_chars_integral<t>, base);
}
//!\endcond

//!\brief Delegates to functions strto[d/f/ld] for floating point value extraction.
template <seqan3::FloatingPoint value_type>
inline std::from_chars_result from_chars_floating_point(char const * first,
                                                        char const * last,
                                                        value_type & value,
                                                        std::chars_format fmt = std::chars_format::general) noexcept
{
    // The locale issue:
    // std::from_chars is documented to be locale independent. The accepted patterns
    // are identical to the one used by strtod in the defailt ("C") locale.
    //
    // The functions strto[d/f/ld] used here are locale dependent but
    // setting the locale manually by std::setlocale is not thread safe.
    // So for the time being this workaround is locale dependent.
    if (*first == '+') // + is permitted in function strto[d/f/ld] but not in from_chars
        return {last, std::errc::invalid_argument};

    float tmp{};
    ptrdiff_t constexpr buffer_size = 100;
    char buffer[buffer_size];

    if (fmt != std::chars_format::general)
    {
        bool exponent_is_present{false};
        for (auto it = first; it != last; ++it)
        {
            if (seqan3::is_char<'e'>(*it) || seqan3::is_char<'E'>(*it))
            {
                exponent_is_present = true;
                break;
            }
        }

        if (fmt == std::chars_format::scientific &&
            !exponent_is_present)
            return {last, std::errc::invalid_argument};

        if (fmt == std::chars_format::fixed      &&
            exponent_is_present)
            return {last, std::errc::invalid_argument};
    }


    // In contrast to std::from_chars, std::strto[f/d/ld] does not treat the second
    // parameter (str_end) as "end of the sequence to parse" but merely as an out
    // parameter to indicate where the parsing ended. Therefore, if [last] does
    // not point to the end of a null-terminated string, a buffer is needed to
    // represent the truncated sequence and ensure correct from_chars functionality.
    char * start;

    if ((*last != '\0' ) || fmt == std::chars_format::hex)
    {
        // If hex format is explicitly expected, the 0x prefix is not allowed in the
        // the original sequence according to the std::from_chars cppreference
        // documentation.
        // In order to use strto[f/d/ld], the prefix must be prepended to achieve
        // correct parsing. This will also automatically lead to an error if the
        // original sequence did contain a 0x prefix and thus reflect the correct
        // requirements of std::from_chars.
        ptrdiff_t offset{0};
        if (fmt == std::chars_format::hex)
        {
            buffer[0] = '0';
            buffer[1] = 'x';
            offset = 2;
        }

        std::copy(first, last, &buffer[offset]);
        buffer[std::min<ptrdiff_t>(buffer_size - offset, last - first)] = '\0';

        start = &buffer[0];
    }
    else
    {
        start = const_cast<char *>(first);
    }

    char * end;

    if constexpr (std::Same<std::remove_reference_t<value_type>, float>)
    {
        tmp = strtof(start, &end);
    }
    if constexpr (std::Same<std::remove_reference_t<value_type>, double>)
    {
        tmp = strtod(start, &end);
    }
    if constexpr (std::Same<std::remove_reference_t<value_type>, long double>)
    {
        tmp = strtold(start, &end);
    }

    last = first + (end - start);

    if (errno == ERANGE)
    {
        return {last, std::errc::result_out_of_range};
    }
    else if (tmp == 0 && end == start)
    {
        return {last, std::errc::invalid_argument};
    }

    // Success.
    value = tmp;
    return {last, {}};
}

} // namespace seqan3::detail
