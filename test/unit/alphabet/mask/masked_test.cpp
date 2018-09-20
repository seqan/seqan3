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
// ==========================================================================

/*!\file
 * \author Joshua Kim <joshua.kim AT fu-berlin.de>
 * \brief Test files for the mask composition.
 */

 #include <gtest/gtest.h>

 #include <seqan3/alphabet/concept_pre.hpp>
 #include <seqan3/alphabet/mask/all.hpp>
 #include <seqan3/alphabet/nucleotide/dna4.hpp>
 #include <seqan3/alphabet/aminoacid/aa20.hpp>

 using namespace seqan3;

 /************** TUPLE INHERITANCE **********************/

 // default/zero construction
 TEST(masked, ctr)
 {
    // Test on dna4.
    [[maybe_unused]] masked<dna4> t1;

    // Test on aa20.
    [[maybe_unused]] masked<aa20> t2;
 }

 // aggregate initialization
 TEST(masked, aggr)
 {
    // Test on dna4.
    [[maybe_unused]]  masked<dna4> t1{dna4::C, mask::MASKED};

    // Test on aa20.
    [[maybe_unused]]  masked<aa20> t2{aa20::W, mask::MASKED};
 }

 // zero initialization
 TEST(masked, zro)
 {
    // Test on dna4.
    masked<dna4> t1{dna4::A, mask::UNMASKED};
    masked<dna4> t2{};
    EXPECT_EQ(t1, t2);

    // Test on aa20.
    masked<aa20> t3{aa20::A, mask::UNMASKED};
    masked<aa20> t4{};
    EXPECT_EQ(t3, t4);
 }

 // copy construction
 TEST(masked, cp_ctr)
 {
    // Test on dna4.
    masked<dna4> t1{dna4::C, mask::MASKED};
    masked<dna4> t2{t1};
    masked<dna4> t3(t1);
    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);

    // Test on aa20.
    masked<aa20> t4{aa20::W, mask::MASKED};
    masked<aa20> t5{t4};
    masked<aa20> t6(t4);
    EXPECT_EQ(t4, t5);
    EXPECT_EQ(t5, t6);
 }

 // move construction
 TEST(masked, mv_ctr)
 {
    // Test on dna4.
    masked<dna4> t0{dna4::C, mask::MASKED};
    masked<dna4> t1{dna4::C, mask::MASKED};
    masked<dna4> t2{std::move(t1)};
    EXPECT_EQ(t2, t0);
    masked<dna4> t3(std::move(t2));
    EXPECT_EQ(t3, t0);

    // Test on aa20.
    masked<aa20> t4{aa20::W, mask::MASKED};
    masked<aa20> t5{aa20::W, mask::MASKED};
    masked<aa20> t6{std::move(t5)};
    EXPECT_EQ(t6, t4);
    masked<aa20> t7(std::move(t6));
    EXPECT_EQ(t7, t4);
 }

 // copy assignment
 TEST(masked, cp_assgn)
 {
    // Test on dna4.
    masked<dna4> t1{dna4::C, mask::MASKED};
    masked<dna4> t2;
    masked<dna4> t3;

    t2 = t1;
    t3 = t1;
    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);

    // Test on aa20.
    masked<aa20> t4{aa20::W, mask::MASKED};
    masked<aa20> t5;
    masked<aa20> t6;

    t5 = t4;
    t6 = t4;
    EXPECT_EQ(t4, t5);
    EXPECT_EQ(t5, t6);
 }

 // move assignment
 TEST(masked, mv_assgn)
 {
    // Test on dna4.
    masked<dna4> t0{dna4::C, mask::MASKED};
    masked<dna4> t1{dna4::C, mask::MASKED};
    masked<dna4> t2;
    masked<dna4> t3;
    t2 = std::move(t1);
    EXPECT_EQ(t2, t0);
    t3 = std::move(t2);
    EXPECT_EQ(t3, t0);

    // Test on aa20.
    masked<aa20> t4{aa20::W, mask::MASKED};
    masked<aa20> t5{aa20::W, mask::MASKED};
    masked<aa20> t6;
    masked<aa20> t7;
    t6 = std::move(t5);
    EXPECT_EQ(t6, t4);
    t7 = std::move(t6);
    EXPECT_EQ(t7, t4);
 }

 // swap
 TEST(masked, swap)
 {
    // Test on dna4.
    masked<dna4> t0{dna4::C, mask::MASKED};
    masked<dna4> t1{dna4::C, mask::MASKED};
    masked<dna4> t2{};
    masked<dna4> t3{};

    std::swap(t1, t2);
    EXPECT_EQ(t2, t0);
    EXPECT_EQ(t1, t3);

    // Test on aa20.
    masked<aa20> t4{aa20::W, mask::MASKED};
    masked<aa20> t5{aa20::W, mask::MASKED};
    masked<aa20> t6{};
    masked<aa20> t7{};

    std::swap(t5, t6);
    EXPECT_EQ(t6, t4);
    EXPECT_EQ(t5, t7);
 }

 // get<1>
 TEST(masked, get_i)
 {
    // Test on dna4.
    masked<dna4> t0{dna4::C, mask::MASKED};

    static_assert(std::is_same_v<decltype(seqan3::get<0>(t0)), dna4 &>);

    EXPECT_EQ(seqan3::get<0>(t0), dna4::C);
    EXPECT_EQ(seqan3::get<1>(t0), mask::MASKED);

    // Test on aa20.
    masked<aa20> t1{aa20::W, mask::MASKED};

    static_assert(std::is_same_v<decltype(seqan3::get<0>(t1)), aa20 &>);

    EXPECT_EQ(seqan3::get<0>(t1), aa20::W);
    EXPECT_EQ(seqan3::get<1>(t1), mask::MASKED);
 }

 // std::get<1>
 TEST(masked, stdget_i)
 {
    // Test on dna4.
    masked<dna4> t0{dna4::C, mask::MASKED};

    static_assert(std::is_same_v<decltype(std::get<0>(t0)), dna4 &>);

    EXPECT_EQ(std::get<0>(t0), dna4::C);
    EXPECT_EQ(std::get<1>(t0), mask::MASKED);

    // Test on aa20.
    masked<aa20> t1{aa20::W, mask::MASKED};

    static_assert(std::is_same_v<decltype(std::get<0>(t1)), aa20 &>);

    EXPECT_EQ(std::get<0>(t1), aa20::W);
    EXPECT_EQ(std::get<1>(t1), mask::MASKED);
 }

 // structured bindings
 TEST(masked, struct_binding)
 {
    // Test on dna4.
    masked<dna4> t0{dna4::C, mask::MASKED};
    auto [ i, l ] = t0;
    static_assert(std::is_same_v<decltype(i), dna4>);
    static_assert(std::is_same_v<decltype(l), mask>);

    EXPECT_EQ(i, dna4::C);
    // EXPECT_EQ(l, mask::MASKED);

    // Test on aa20.
    masked<aa20> t1{aa20::W, mask::MASKED};
    auto [ j, k ] = t1;

    static_assert(std::is_same_v<decltype(j), aa20>);
    static_assert(std::is_same_v<decltype(k), mask>);

    EXPECT_EQ(j, aa20::W);
    // EXPECT_EQ(k, mask::MASKED);
 }

 // get<type>
 TEST(masked, get_type)
 {
    // Test on dna4.
    masked<dna4> t0{dna4::C, mask::MASKED};

    EXPECT_EQ(seqan3::get<dna4>(t0), dna4::C);

    // Test on aa20.
    masked<aa20> t1{aa20::W, mask::MASKED};

    EXPECT_EQ(seqan3::get<aa20>(t1), aa20::W);
 }

 // std::get<type>
 TEST(masked, stdget_type)
 {
    // Test on dna4.
    masked<dna4> t0{dna4::C, mask::MASKED};

    EXPECT_EQ(std::get<dna4>(t0), dna4::C);

    // Test on aa20.
    masked<aa20> t1{aa20::W, mask::MASKED};

    EXPECT_EQ(std::get<aa20>(t1), aa20::W);
 }

 // std::tuple_element
 TEST(masked, tuple_element)
 {
    // Test on dna4.
    using pt = masked<dna4>;

    static_assert(std::is_same_v<std::tuple_element_t<0, pt>, dna4>);
    static_assert(std::tuple_size_v<pt> == 2);

    // Test on aa20.
    using aapt = masked<aa20>;

    static_assert(std::is_same_v<std::tuple_element_t<0, aapt>, aa20>);
    static_assert(std::tuple_size_v<aapt> == 2);
 }

 // type deduction
 TEST(masked, type_deduce)
 {
    // Test on dna4.
    masked t0{dna4::C, mask::MASKED};
    using pt = decltype(t0);

    static_assert(std::is_same_v<std::tuple_element_t<0, pt>, dna4>);
    static_assert(std::tuple_size_v<pt> == 2);

    // Test on aa20.
    masked t1{aa20::W, mask::MASKED};
    using aapt = decltype(t1);

    static_assert(std::is_same_v<std::tuple_element_t<0, aapt>, aa20>);
    static_assert(std::tuple_size_v<aapt> == 2);
 }

 // explicit cast to element
 TEST(masked, cast_to_element)
 {
    // Test on dna4.
    masked<dna4> t0{dna4::C, mask::MASKED};

    auto d = static_cast<dna4>(t0);
    auto q = static_cast<mask>(t0);
    static_assert(std::is_same_v<decltype(d), dna4>);

    EXPECT_EQ(d, dna4::C);
    EXPECT_EQ(q, mask::MASKED);

    // Test on aa20.
    masked<aa20> t1{aa20::W, mask::MASKED};

    auto e = static_cast<aa20>(t1);
    auto r = static_cast<mask>(t1);
    static_assert(std::is_same_v<decltype(e), aa20>);

    EXPECT_EQ(e, aa20::W);
    EXPECT_EQ(r, mask::MASKED);
 }

 // comparison operators
 TEST(masked, cmp)
 {
    // Test on dna4.
    masked<dna4> t0{dna4::C, mask::UNMASKED};
    masked<dna4> t1{dna4::C, mask::MASKED};
    masked<dna4> t2{dna4::G, mask::MASKED};

    EXPECT_LT(t0, t1);
    EXPECT_LE(t0, t1);
    EXPECT_LE(t1, t1);
    EXPECT_EQ(t1, t1);
    EXPECT_GE(t1, t1);
    EXPECT_GE(t2, t1);
    EXPECT_GT(t2, t1);

    // Test on aa20.
    masked<aa20> t3{aa20::W, mask::UNMASKED};
    masked<aa20> t4{aa20::W, mask::MASKED};
    masked<aa20> t5{aa20::Y, mask::MASKED};

    EXPECT_LT(t3, t4);
    EXPECT_LE(t3, t4);
    EXPECT_LE(t4, t4);
    EXPECT_EQ(t4, t4);
    EXPECT_GE(t4, t4);
    EXPECT_GE(t5, t4);
    EXPECT_GT(t5, t4);
 }

/************** SEMI ALPHABET and ALPHABET concept **********************/
TEST(masked, rank_type)
{
    // Test on dna4.
    EXPECT_TRUE((std::is_same_v<underlying_rank_t<masked<dna4>>,
                               uint8_t>));
    // Test on aa20.
    EXPECT_TRUE((std::is_same_v<underlying_rank_t<masked<aa20>>,
                               uint8_t>));
}

TEST(masked, char_type)
{
    // Test on dna4.
    EXPECT_TRUE((std::is_same_v<underlying_char_t<masked<dna4>>,
                               underlying_char_t<dna4>>));

    // Test on aa20.
    EXPECT_TRUE((std::is_same_v<underlying_char_t<masked<aa20>>,
                               underlying_char_t<aa20>>));
}

TEST(masked, alphabet_size_v)
{
    // Test on dna4.
    EXPECT_EQ((alphabet_size_v<masked<dna4>>),
              (alphabet_size_v<dna4> * alphabet_size_v<mask>));

    // Test on aa20.
    EXPECT_EQ((alphabet_size_v<masked<aa20>>),
              (alphabet_size_v<aa20> * alphabet_size_v<mask>));
}

TEST(masked, to_rank)
{
    // Test on dna4.
    masked<dna4> t0{dna4::C, mask::MASKED};
    EXPECT_EQ(to_rank(std::get<0>(t0)), 1);
    // EXPECT_EQ(to_rank(std::get<1>(t0)), 1);
    // EXPECT_EQ(to_rank(t0),
    //           to_rank(std::get<0>(t0)) +
    //           alphabet_size_v<dna4> * to_rank(std::get<1>(t0)));

    // Test on aa20.
    masked<aa20> t1{aa20::A, mask::UNMASKED};
    EXPECT_EQ(to_rank(std::get<0>(t1)), 0);
    // EXPECT_EQ(to_rank(std::get<1>(t1)), 0);
    // EXPECT_EQ(to_rank(t1),
    //           to_rank(std::get<0>(t1)) +
    //           alphabet_size_v<aa20> * to_rank(std::get<1>(t1)));
}

TEST(masked, assign_rank)
{
    // Test on dna4.
    using type = masked<dna4>;

    type t0{};

    for (underlying_rank_t<type> i = 0; i < alphabet_size_v<type>; ++i)
    {
        assign_rank(t0, i);
        EXPECT_EQ(to_rank(t0), i);
    }

    // Test on aa20.
    using aatype = masked<aa20>;

    aatype t1{};

    for (underlying_rank_t<aatype> j = 0; j < alphabet_size_v<aatype>; ++j)
    {
        assign_rank(t1, j);
        EXPECT_EQ(to_rank(t1), j);
    }
}

TEST(masked, to_char)
{
    // Test on dna4.
    masked<dna4> t0{dna4::C, mask::UNMASKED};
    EXPECT_EQ(to_char(std::get<0>(t0)), 'C');
    EXPECT_EQ(to_char(t0), 'C');
    t0 = mask::MASKED;
    EXPECT_EQ(to_char(std::get<0>(t0)), 'C');
    EXPECT_EQ(to_char(t0), 'c');

    // Test on aa20.
    masked<aa20> t1{aa20::W, mask::UNMASKED};
    EXPECT_EQ(to_char(std::get<0>(t1)), 'W');
    EXPECT_EQ(to_char(t1), 'W');
    t1 = mask::MASKED;
    EXPECT_EQ(to_char(std::get<0>(t1)), 'W');
    EXPECT_EQ(to_char(t1), 'w');
}

TEST(masked, assign_char)
{
    // Test on dna4.
    using type = masked<dna4>;
    type t0{dna4::C, mask::MASKED};
    assign_char(t0, 'A');
    EXPECT_EQ(to_char(t0), 'A');
    assign_char(t0, 'C');
    EXPECT_EQ(to_char(t0), 'C');
    assign_char(t0, 'G');
    EXPECT_EQ(to_char(t0), 'G');
    assign_char(t0, 'T');
    EXPECT_EQ(to_char(t0), 'T');
    assign_char(t0, 'N');
    EXPECT_EQ(to_char(t0), 'A');
    assign_char(t0, 'a');
    EXPECT_EQ(to_char(t0), 'a');
    assign_char(t0, 'c');
    EXPECT_EQ(to_char(t0), 'c');
    assign_char(t0, 'g');
    EXPECT_EQ(to_char(t0), 'g');
    assign_char(t0, 't');
    EXPECT_EQ(to_char(t0), 't');

    // Test on aa20.
    using aatype = masked<aa20>;
    aatype t1{aa20::W, mask::MASKED};
    assign_char(t1, 'A');
    EXPECT_EQ(to_char(t1), 'A');
    assign_char(t1, 'C');
    EXPECT_EQ(to_char(t1), 'C');
    assign_char(t1, 'D');
    EXPECT_EQ(to_char(t1), 'D');
    assign_char(t1, 'E');
    EXPECT_EQ(to_char(t1), 'E');
    assign_char(t1, 'F');
    EXPECT_EQ(to_char(t1), 'F');
    assign_char(t1, 'G');
    EXPECT_EQ(to_char(t1), 'G');
    assign_char(t1, 'H');
    EXPECT_EQ(to_char(t1), 'H');
    assign_char(t1, 'I');
    EXPECT_EQ(to_char(t1), 'I');
    assign_char(t1, 'K');
    EXPECT_EQ(to_char(t1), 'K');
    assign_char(t1, 'L');
    EXPECT_EQ(to_char(t1), 'L');
    assign_char(t1, 'M');
    EXPECT_EQ(to_char(t1), 'M');
    assign_char(t1, 'N');
    EXPECT_EQ(to_char(t1), 'N');
    assign_char(t1, 'P');
    EXPECT_EQ(to_char(t1), 'P');
    assign_char(t1, 'Q');
    EXPECT_EQ(to_char(t1), 'Q');
    assign_char(t1, 'R');
    EXPECT_EQ(to_char(t1), 'R');
    assign_char(t1, 'S');
    EXPECT_EQ(to_char(t1), 'S');
    assign_char(t1, 'T');
    EXPECT_EQ(to_char(t1), 'T');
    assign_char(t1, 'V');
    EXPECT_EQ(to_char(t1), 'V');
    assign_char(t1, 'W');
    EXPECT_EQ(to_char(t1), 'W');
    assign_char(t1, 'Y');
    EXPECT_EQ(to_char(t1), 'Y');
    assign_char(t1, 'B');
    EXPECT_EQ(to_char(t1), 'D');
    assign_char(t1, 'J');
    EXPECT_EQ(to_char(t1), 'L');
    assign_char(t1, 'O');
    EXPECT_EQ(to_char(t1), 'L');
    assign_char(t1, 'U');
    EXPECT_EQ(to_char(t1), 'C');
    assign_char(t1, 'X');
    EXPECT_EQ(to_char(t1), 'S');
    assign_char(t1, 'Z');
    EXPECT_EQ(to_char(t1), 'E');
    assign_char(t1, to_char(aa20::TERMINATOR));
    EXPECT_EQ(to_char(t1), 'W');
    assign_char(t1, to_char(aa20::UNKNOWN));
    EXPECT_EQ(to_char(t1), 'S');
    assign_char(t1, 'a');
    EXPECT_EQ(to_char(t1), 'a');
    assign_char(t1, 'c');
    EXPECT_EQ(to_char(t1), 'c');
    assign_char(t1, 'd');
    EXPECT_EQ(to_char(t1), 'd');
    assign_char(t1, 'e');
    EXPECT_EQ(to_char(t1), 'e');
    assign_char(t1, 'f');
    EXPECT_EQ(to_char(t1), 'f');
    assign_char(t1, 'g');
    EXPECT_EQ(to_char(t1), 'g');
    assign_char(t1, 'h');
    EXPECT_EQ(to_char(t1), 'h');
    assign_char(t1, 'i');
    EXPECT_EQ(to_char(t1), 'i');
    assign_char(t1, 'k');
    EXPECT_EQ(to_char(t1), 'k');
    assign_char(t1, 'l');
    EXPECT_EQ(to_char(t1), 'l');
    assign_char(t1, 'm');
    EXPECT_EQ(to_char(t1), 'm');
    assign_char(t1, 'n');
    EXPECT_EQ(to_char(t1), 'n');
    assign_char(t1, 'p');
    EXPECT_EQ(to_char(t1), 'p');
    assign_char(t1, 'q');
    EXPECT_EQ(to_char(t1), 'q');
    assign_char(t1, 'r');
    EXPECT_EQ(to_char(t1), 'r');
    assign_char(t1, 's');
    EXPECT_EQ(to_char(t1), 's');
    assign_char(t1, 't');
    EXPECT_EQ(to_char(t1), 't');
    assign_char(t1, 'v');
    EXPECT_EQ(to_char(t1), 'v');
    assign_char(t1, 'w');
    EXPECT_EQ(to_char(t1), 'w');
    assign_char(t1, 'y');
    EXPECT_EQ(to_char(t1), 'y');
    assign_char(t1, 'b');
    EXPECT_EQ(to_char(t1), 'd');
    assign_char(t1, 'j');
    EXPECT_EQ(to_char(t1), 'l');
    assign_char(t1, 'o');
    EXPECT_EQ(to_char(t1), 'l');
    assign_char(t1, 'u');
    EXPECT_EQ(to_char(t1), 'c');
    assign_char(t1, 'x');
    EXPECT_EQ(to_char(t1), 's');
    assign_char(t1, 'z');
    EXPECT_EQ(to_char(t1), 'e');
    assign_char(t1, tolower(to_char(aa20::TERMINATOR)));
    EXPECT_EQ(to_char(t1), 'w');
    assign_char(t1, tolower(to_char(aa20::UNKNOWN)));
    EXPECT_EQ(to_char(t1), 's');
}
