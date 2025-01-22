// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Sara Hetzel <sara.hetzel AT fu-berlin.de>
 * \brief Provides translation details for nucleotide to aminoacid translation.
 */

#pragma once

#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/aminoacid/translation_genetic_code.hpp>
#include <seqan3/alphabet/nucleotide/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna15.hpp>

// clang-format off
namespace seqan3::detail
{
/*!\brief Generic translation table for canonical genetic code.
 * \tparam nucl_type The type of input nucleotides.
 *
 * \details  All nucleotides are casted to dna15 and the dna15 translation table is used for translation.
 */
template <typename nucl_type, seqan3::genetic_code gc = seqan3::genetic_code::canonical, typename void_type = void>
struct translation_table
{
    //!\brief Holds the generic translation table.
    static constexpr std::array<std::array<std::array<aa27, alphabet_size<nucl_type>>, alphabet_size<nucl_type>>,
                                                      alphabet_size<nucl_type>> value
    {
        [] () constexpr
        {
            std::array<std::array<std::array<aa27, alphabet_size<nucl_type>>,
                                             alphabet_size<nucl_type>>, alphabet_size<nucl_type>> table{};

            using size_t = std::remove_const_t<decltype(alphabet_size<nucl_type>)>;
            for (size_t i = 0; i < alphabet_size<nucl_type>; ++i)
            {
                dna15 n1(assign_rank_to(i, nucl_type{}));
                for (size_t j = 0; j < alphabet_size<nucl_type>; ++j)
                {
                    dna15 n2(assign_rank_to(j, nucl_type{}));
                    for (size_t k = 0; k < alphabet_size<nucl_type>; ++k)
                    {
                        dna15 n3(assign_rank_to(k, nucl_type{}));
                        table[i][j][k] = translation_table<dna15, gc, void_type>::value[to_rank(n1)][to_rank(n2)][to_rank(n3)];
                    }
                }
            }
            return table;
        } ()
    };
};

//!\brief Translation table for canonical genetic code and dna15 alphabet.
template <typename void_type>
struct translation_table<dna15, seqan3::genetic_code::canonical, void_type>
{
    //!\brief Holds the translation table for canonical genetic code and nucl16 alphabet.
    static constexpr aa27 value[dna15::alphabet_size][dna15::alphabet_size][dna15::alphabet_size]
    {
        { // a??
            // a,        b,        c,        d,        g,        h,        k,        m,        n,        r,        s,        t,        v,        w,        y
            { 'K'_aa27, 'X'_aa27, 'N'_aa27, 'X'_aa27, 'K'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'K'_aa27, 'X'_aa27, 'N'_aa27, 'X'_aa27, 'X'_aa27, 'N'_aa27 }, // aa?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // ab?
            { 'T'_aa27, 'T'_aa27, 'T'_aa27, 'T'_aa27, 'T'_aa27, 'T'_aa27, 'T'_aa27, 'T'_aa27, 'T'_aa27, 'T'_aa27, 'T'_aa27, 'T'_aa27, 'T'_aa27, 'T'_aa27, 'T'_aa27 }, // ac?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // ad?
            { 'R'_aa27, 'X'_aa27, 'S'_aa27, 'X'_aa27, 'R'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'R'_aa27, 'X'_aa27, 'S'_aa27, 'X'_aa27, 'X'_aa27, 'S'_aa27 }, // ag?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // ah?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // ak?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // am?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // an?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // ar?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // as?
            { 'I'_aa27, 'X'_aa27, 'I'_aa27, 'X'_aa27, 'M'_aa27, 'I'_aa27, 'X'_aa27, 'I'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'I'_aa27, 'X'_aa27, 'I'_aa27, 'I'_aa27 }, // at?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // av?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // aw?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }  // ay?
    }, { // b??
            // a,        b,        c,        d,        g,        h,        k,        m,        n,        r,        s,        t,        v,        w,        y
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // ba?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // bb?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // bc?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // bd?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // bg?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // bh?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // bk?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // bm?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // bn?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // br?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // bs?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // bt?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // bv?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // bw?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }  // by?
        }, { // c??
            // a,        b,        c,        d,        g,        h,        k,        m,        n,        r,        s,        t,        v,        w,        y
            { 'Q'_aa27, 'X'_aa27, 'H'_aa27, 'X'_aa27, 'Q'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'Q'_aa27, 'X'_aa27, 'H'_aa27, 'X'_aa27, 'X'_aa27, 'H'_aa27 }, // ca?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // cb?
            { 'P'_aa27, 'P'_aa27, 'P'_aa27, 'P'_aa27, 'P'_aa27, 'P'_aa27, 'P'_aa27, 'P'_aa27, 'P'_aa27, 'P'_aa27, 'P'_aa27, 'P'_aa27, 'P'_aa27, 'P'_aa27, 'P'_aa27 }, // cc?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // cd?
            { 'R'_aa27, 'R'_aa27, 'R'_aa27, 'R'_aa27, 'R'_aa27, 'R'_aa27, 'R'_aa27, 'R'_aa27, 'R'_aa27, 'R'_aa27, 'R'_aa27, 'R'_aa27, 'R'_aa27, 'R'_aa27, 'R'_aa27 }, // cg?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // ch?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // ck?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // cm?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // cn?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // cr?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // cs?
            { 'L'_aa27, 'L'_aa27, 'L'_aa27, 'L'_aa27, 'L'_aa27, 'L'_aa27, 'L'_aa27, 'L'_aa27, 'L'_aa27, 'L'_aa27, 'L'_aa27, 'L'_aa27, 'L'_aa27, 'L'_aa27, 'L'_aa27 }, // ct?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // cv?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // cw?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }  // cy?
        }, { // d??
            // a,        b,        c,        d,        g,        h,        k,        m,        n,        r,        s,        t,        v,        w,        y
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // da?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // db?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // dc?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // dd?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // dg?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // dh?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // dk?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // dm?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // dn?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // dr?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // ds?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // dt?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // dv?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // dw?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }  // dy?
        }, { // g??
            // a,        b,        c,        d,        g,        h,        k,        m,        n,        r,        s,        t,        v,        w,        y
            { 'E'_aa27, 'X'_aa27, 'D'_aa27, 'X'_aa27, 'E'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'E'_aa27, 'X'_aa27, 'D'_aa27, 'X'_aa27, 'X'_aa27, 'D'_aa27 }, // ga?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // gb?
            { 'A'_aa27, 'A'_aa27, 'A'_aa27, 'A'_aa27, 'A'_aa27, 'A'_aa27, 'A'_aa27, 'A'_aa27, 'A'_aa27, 'A'_aa27, 'A'_aa27, 'A'_aa27, 'A'_aa27, 'A'_aa27, 'A'_aa27 }, // gc?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // gd?
            { 'G'_aa27, 'G'_aa27, 'G'_aa27, 'G'_aa27, 'G'_aa27, 'G'_aa27, 'G'_aa27, 'G'_aa27, 'G'_aa27, 'G'_aa27, 'G'_aa27, 'G'_aa27, 'G'_aa27, 'G'_aa27, 'G'_aa27 }, // gg?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // gh?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // gk?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // gm?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // gn?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // gr?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // gs?
            { 'V'_aa27, 'V'_aa27, 'V'_aa27, 'V'_aa27, 'V'_aa27, 'V'_aa27, 'V'_aa27, 'V'_aa27, 'V'_aa27, 'V'_aa27, 'V'_aa27, 'V'_aa27, 'V'_aa27, 'V'_aa27, 'V'_aa27 }, // gt?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // gv?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // gw?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }  // gy?
        }, { // h??
            // a,        b,        c,        d,        g,        h,        k,        m,        n,        r,        s,        t,        v,        w,        y
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // ha?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // hb?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // hc?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // hd?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // hg?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // hh?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // hk?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // hm?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // hn?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // hr?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // hs?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // ht?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // hv?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // hw?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }  // hy?
        }, { // k??
            // a,        b,        c,        d,        g,        h,        k,        m,        n,        r,        s,        t,        v,        w,        y
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // ka?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // kb?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // kc?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // kd?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // kg?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // kh?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // kk?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // km?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // kn?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // kr?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // ks?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // kt?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // kv?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // kw?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }  // ky?
        }, { // m??
            // a,        b,        c,        d,        g,        h,        k,        m,        n,        r,        s,        t,        v,        w,        y
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // ma?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // mb?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // mc?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // md?
            { 'R'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'R'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'R'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // mg?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // mh?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // mk?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // mm?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // mn?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // mr?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // ms?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // mt?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // mv?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // mw?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }  // my?
        }, { // n??
            // a,        b,        c,        d,        g,        h,        k,        m,        n,        r,        s,        t,        v,        w,        y
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // na?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // nb?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // nc?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // nd?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // ng?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // nh?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // nk?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // nm?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // nn?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // nr?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // ns?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // nt?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // nv?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // nw?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }  // ny?
        }, { // r??
            // a,        b,        c,        d,        g,        h,        k,        m,        n,        r,        s,        t,        v,        w,        y
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // ra?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // rb?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // rc?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // rd?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // rg?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // rh?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // rk?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // rm?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // rn?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // rr?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // rs?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // rt?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // rv?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // rw?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }  // ry?
        }, { // s??
            // a,        b,        c,        d,        g,        h,        k,        m,        n,        r,        s,        t,        v,        w,        y
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // sa?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // sb?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // sc?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // sd?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // sg?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // sh?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // sk?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // sm?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // sn?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // sr?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // ss?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // st?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // sv?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // sw?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }  // sy?
        }, { // t??
            // a,        b,        c,        d,        g,        h,        k,        m,        n,        r,        s,        t,        v,        w,        y
            { '*'_aa27, 'X'_aa27, 'Y'_aa27, 'X'_aa27, '*'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, '*'_aa27, 'X'_aa27, 'Y'_aa27, 'X'_aa27, 'X'_aa27, 'Y'_aa27 }, // ta?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // tb?
            { 'S'_aa27, 'S'_aa27, 'S'_aa27, 'S'_aa27, 'S'_aa27, 'S'_aa27, 'S'_aa27, 'S'_aa27, 'S'_aa27, 'S'_aa27, 'S'_aa27, 'S'_aa27, 'S'_aa27, 'S'_aa27, 'S'_aa27 }, // tc?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // td?
            { '*'_aa27, 'X'_aa27, 'C'_aa27, 'X'_aa27, 'W'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'C'_aa27, 'X'_aa27, 'X'_aa27, 'C'_aa27 }, // tg?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // th?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // tk?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // tm?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // tn?
            { '*'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // tr?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // ts?
            { 'L'_aa27, 'X'_aa27, 'F'_aa27, 'X'_aa27, 'L'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'L'_aa27, 'X'_aa27, 'F'_aa27, 'X'_aa27, 'X'_aa27, 'F'_aa27 }, // tt?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // tv?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // tw?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }  // ty?
        }, { // v??
            // a,        b,        c,        d,        g,        h,        k,        m,        n,        r,        s,        t,        v,        w,        y
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // va?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // vb?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // vc?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // vd?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // vg?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // vh?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // vk?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // vm?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // vn?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // vr?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // vs?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // vt?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // vv?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // vw?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }  // vy?
        }, { // w??
            // a,        b,        c,        d,        g,        h,        k,        m,        n,        r,        s,        t,        v,        w,        y
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // wa?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // wb?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // wc?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // wd?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // wg?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // wh?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // wk?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // wm?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // wn?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // wr?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // ws?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // wt?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // wv?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // ww?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }  // wy?
        }, { // y??
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // ya?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // yb?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // yc?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // yd?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // yg?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // yh?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // yk?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // ym?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // yn?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // yr?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // ys?
            { 'L'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'L'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'L'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // yt?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // yv?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }, // yw?
            { 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27, 'X'_aa27 }  // yy?
        }
    };
};

} // namespace seqan3::detail
// clang-format on
