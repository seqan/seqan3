// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin 
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik 
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Sara Hetzel <sara.hetzel AT fu-berlin.de>
 * \brief Contains translation details for nucleotide to aminoacid translation.
 */

#pragma once

#include <seqan3/alphabet/nucleotide/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/aminoacid/translation_genetic_code.hpp>

namespace seqan3::detail
{
/*!\brief Generic translation table for canonical genetic code.
 * \tparam nucl_type The type of input nucleotides.
 *
 * \details  All nucleotides are casted to dna15 and the dna15 translation table is used for translation.
 */
template <typename nucl_type, seqan3::genetic_code gc = seqan3::genetic_code::CANONICAL, typename void_type = void>
struct translation_table
{
    //!\brief Holds the generic translation table.
    static constexpr std::array<std::array<std::array<aa27, alphabet_size_v<nucl_type>>, alphabet_size_v<nucl_type>>,
                                                      alphabet_size_v<nucl_type>> VALUE
    {
        [] () constexpr
        {
            std::array<std::array<std::array<aa27, alphabet_size_v<nucl_type>>,
                                             alphabet_size_v<nucl_type>>, alphabet_size_v<nucl_type>> table{};

            using size_t = std::remove_const_t<decltype(alphabet_size_v<nucl_type>)>;
            for (size_t i = 0; i < alphabet_size_v<nucl_type>; ++i)
            {
                dna15 n1(assign_rank(nucl_type{}, i));
                for (size_t j = 0; j < alphabet_size_v<nucl_type>; ++j)
                {
                    dna15 n2(assign_rank(nucl_type{}, j));
                    for (size_t k = 0; k < alphabet_size_v<nucl_type>; ++k)
                    {
                        dna15 n3(assign_rank(nucl_type{}, k));
                        table[i][j][k] = translation_table<dna15, gc, void_type>::VALUE[to_rank(n1)][to_rank(n2)][to_rank(n3)];
                    }
                }
            }
            return table;
        } ()
    };
};

//!\brief Translation table for canonical genetic code and dna15 alphabet.
template <typename void_type>
struct translation_table<dna15, seqan3::genetic_code::CANONICAL, void_type>
{
    //!\brief Holds the translation table for canonical genetic code and nucl16 alphabet.
    static constexpr aa27 VALUE[dna15::value_size][dna15::value_size][dna15::value_size]
    {
        { // a??
            // a,               b,       c,       d,       g,                h,       k,       m,       n,       r,                s,       t,       v,       w,       y
            { aa27::K,          aa27::X, aa27::N, aa27::X, aa27::K,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::K,          aa27::X, aa27::N, aa27::X, aa27::X, aa27::N }, // aa?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ab?
            { aa27::T,          aa27::T, aa27::T, aa27::T, aa27::T,          aa27::T, aa27::T, aa27::T, aa27::T, aa27::T,          aa27::T, aa27::T, aa27::T, aa27::T, aa27::T }, // ac?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ad?
            { aa27::R,          aa27::X, aa27::S, aa27::X, aa27::R,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::R,          aa27::X, aa27::S, aa27::X, aa27::X, aa27::S }, // ag?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ah?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ak?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // am?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // an?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ar?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // as?
            { aa27::I,          aa27::X, aa27::I, aa27::X, aa27::M,          aa27::I, aa27::X, aa27::I, aa27::X, aa27::X,          aa27::X, aa27::I, aa27::X, aa27::I, aa27::I }, // at?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // av?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // aw?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }  // ay?
    }, { // b??
            // a,               b,       c,       d,       g,                h,       k,       m,       n,       r,                s,       t,       v,       w,       y
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ba?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // bb?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // bc?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // bd?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // bg?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // bh?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // bk?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // bm?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // bn?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // br?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // bs?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // bt?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // bv?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // bw?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }  // by?
        }, { // c??
            // a,               b,       c,       d,       g,                h,       k,       m,       n,       r,                s,       t,       v,       w,       y
            { aa27::Q,          aa27::X, aa27::H, aa27::X, aa27::Q,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::Q,          aa27::X, aa27::H, aa27::X, aa27::X, aa27::H }, // ca?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // cb?
            { aa27::P,          aa27::P, aa27::P, aa27::P, aa27::P,          aa27::P, aa27::P, aa27::P, aa27::P, aa27::P,          aa27::P, aa27::P, aa27::P, aa27::P, aa27::P }, // cc?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // cd?
            { aa27::R,          aa27::R, aa27::R, aa27::R, aa27::R,          aa27::R, aa27::R, aa27::R, aa27::R, aa27::R,          aa27::R, aa27::R, aa27::R, aa27::R, aa27::R }, // cg?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ch?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ck?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // cm?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // cn?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // cr?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // cs?
            { aa27::L,          aa27::L, aa27::L, aa27::L, aa27::L,          aa27::L, aa27::L, aa27::L, aa27::L, aa27::L,          aa27::L, aa27::L, aa27::L, aa27::L, aa27::L }, // ct?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // cv?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // cw?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }  // cy?
        }, { // d??
            // a,               b,       c,       d,       g,                h,       k,       m,       n,       r,                s,       t,       v,       w,       y
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // da?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // db?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // dc?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // dd?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // dg?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // dh?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // dk?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // dm?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // dn?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // dr?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ds?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // dt?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // dv?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // dw?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }  // dy?

        }, { // g??
            // a,               b,       c,       d,       g,                h,       k,       m,       n,       r,                s,       t,       v,       w,       y
            { aa27::E,          aa27::X, aa27::D, aa27::X, aa27::E,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::E,          aa27::X, aa27::D, aa27::X, aa27::X, aa27::D }, // ga?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // gb?
            { aa27::A,          aa27::A, aa27::A, aa27::A, aa27::A,          aa27::A, aa27::A, aa27::A, aa27::A, aa27::A,          aa27::A, aa27::A, aa27::A, aa27::A, aa27::A }, // gc?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // gd?
            { aa27::G,          aa27::G, aa27::G, aa27::G, aa27::G,          aa27::G, aa27::G, aa27::G, aa27::G, aa27::G,          aa27::G, aa27::G, aa27::G, aa27::G, aa27::G }, // gg?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // gh?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // gk?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // gm?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // gn?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // gr?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // gs?
            { aa27::V,          aa27::V, aa27::V, aa27::V, aa27::V,          aa27::V, aa27::V, aa27::V, aa27::V, aa27::V,          aa27::V, aa27::V, aa27::V, aa27::V, aa27::V }, // gt?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // gv?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // gw?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }  // gy?
        }, { // h??
            // a,               b,       c,       d,       g,                h,       k,       m,       n,       r,                s,       t,       v,       w,       y
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ha?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // hb?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // hc?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // hd?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // hg?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // hh?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // hk?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // hm?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // hn?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // hr?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // hs?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ht?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // hv?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // hw?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }  // hy?

        }, { // k??
            // a,               b,       c,       d,       g,                h,       k,       m,       n,       r,                s,       t,       v,       w,       y
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ka?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // kb?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // kc?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // kd?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // kg?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // kh?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // kk?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // km?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // kn?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // kr?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ks?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // kt?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // kv?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // kw?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }  // ky?

        }, { // m??
            // a,               b,       c,       d,       g,                h,       k,       m,       n,       r,                s,       t,       v,       w,       y
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ma?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // mb?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // mc?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // md?
            { aa27::R,          aa27::X, aa27::X, aa27::X, aa27::R,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::R,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // mg?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // mh?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // mk?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // mm?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // mn?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // mr?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ms?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // mt?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // mv?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // mw?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }  // my?


        }, { // n??
            // a,               b,       c,       d,       g,                h,       k,       m,       n,       r,                s,       t,       v,       w,       y
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // nb?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // nc?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // nd?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ng?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // nh?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // nk?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // nm?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // nn?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // nr?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ns?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // nt?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // nv?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // nw?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }  // ny?
        }, { // r??
            // a,               b,       c,       d,       g,                h,       k,       m,       n,       r,                s,       t,       v,       w,       y
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // rb?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // rc?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // rd?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // rg?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // rh?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // rk?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // rm?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // rn?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // rr?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // rs?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // rt?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // rv?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // rw?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }  // ry?
        }, { // s??
            // a,               b,       c,       d,       g,                h,       k,       m,       n,       r,                s,       t,       v,       w,       y
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // sb?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // sc?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // sd?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // sg?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // sh?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // sk?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // sm?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // sn?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // sr?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ss?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // st?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // sv?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // sw?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }  // sy?
        }, { // t??
            // a,               b,       c,       d,       g,                h,       k,       m,       n,       r,                s,       t,       v,       w,       y
            { aa27::TERMINATOR, aa27::X, aa27::Y, aa27::X, aa27::TERMINATOR, aa27::X, aa27::X, aa27::X, aa27::X, aa27::TERMINATOR, aa27::X, aa27::Y, aa27::X, aa27::X, aa27::Y }, // ta?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // tb?
            { aa27::S,          aa27::S, aa27::S, aa27::S, aa27::S,          aa27::S, aa27::S, aa27::S, aa27::S, aa27::S,          aa27::S, aa27::S, aa27::S, aa27::S, aa27::S }, // tc?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // td?
            { aa27::TERMINATOR, aa27::X, aa27::C, aa27::X, aa27::W,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::C, aa27::X, aa27::C }, // tg?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // th?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // tk?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // tm?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // tn?
            { aa27::TERMINATOR, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // tr?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ts?
            { aa27::L,          aa27::X, aa27::F, aa27::X, aa27::X,          aa27::L, aa27::X, aa27::X, aa27::X, aa27::L,          aa27::X, aa27::F, aa27::X, aa27::X, aa27::F }, // tt?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // tv?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // tw?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }  // ty?
        }, { // v??
            // a,               b,       c,       d,       g,                h,       k,       m,       n,       r,                s,       t,       v,       w,       y
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // va?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // vb?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // vc?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // vd?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // vg?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // vh?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // vk?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // vm?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // vn?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // vr?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // vs?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // vt?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // vv?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // vw?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }  // vy?
        }, { // w??
            // a,               b,       c,       d,       g,                h,       k,       m,       n,       r,                s,       t,       v,       w,       y
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // wa?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // wb?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // wc?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // wd?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // wg?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // wh?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // wk?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // wm?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // wn?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // wr?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ws?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // wt?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // wv?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ww?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }  // wy?
        }, { // y??
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ya?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // yb?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // yc?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // yd?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // yg?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // yh?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // yk?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ym?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // yn?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // yr?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // ys?
            { aa27::L,          aa27::X, aa27::X, aa27::X, aa27::L,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::L,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // yt?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // yv?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }, // yw?
            { aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X,          aa27::X, aa27::X, aa27::X, aa27::X, aa27::X }  // yy?
        }
    };
};

}
