// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <seqan3/alignment/pairwise/align_result_selector.hpp>
#include <seqan3/alignment/pairwise/alignment_result.hpp>
#include <seqan3/alignment/pairwise/edit_distance_unbanded.hpp>
#include <seqan3/test/seqan2.hpp>

#ifdef SEQAN3_HAS_SEQAN2
#    include <seqan/find.h>
#endif

namespace seqan3::test
{

// ----------------------------------------------------------------------------
//  edit distance pairwise alignment (seqan3)
// ----------------------------------------------------------------------------

struct edit_distance_algorithm
{
    template <typename sequence1_t, typename sequence2_t, typename edit_distance_cfg_t>
    using alignment_result_type =
        seqan3::alignment_result<typename seqan3::detail::align_result_selector<sequence1_t const &,
                                                                                sequence2_t const &,
                                                                                edit_distance_cfg_t>::type>;

    template <typename sequence1_t, typename sequence2_t, typename edit_distance_cfg_t>
    static auto select(edit_distance_cfg_t edit_distance_cfg)
    {
        using seqan3::get;
        auto method_global = get<seqan3::align_cfg::method_global>(edit_distance_cfg);
        assert(method_global.free_end_gaps_sequence2_leading == false);
        assert(method_global.free_end_gaps_sequence2_trailing == false);
        assert(method_global.free_end_gaps_sequence1_leading == method_global.free_end_gaps_sequence1_trailing);

        // NOTE: returning lambdas resulted in faster benchmark times than simply returning
        // `algorithm_impl<std::{false,true}_type>::execute` as function pointer of type
        // `alignment_result_t (*) (sequence1_t const &, sequence2_t const &, edit_distance_cfg_t)`
        return method_global.free_end_gaps_sequence1_leading ?
               [] (sequence1_t const & sequence1, sequence2_t const & sequence2, edit_distance_cfg_t edit_distance_cfg)
               {
                   return algorithm_impl<std::true_type>::execute(sequence1, sequence2, edit_distance_cfg);
               } :
               [] (sequence1_t const & sequence1, sequence2_t const & sequence2, edit_distance_cfg_t edit_distance_cfg)
               {
                   return algorithm_impl<std::false_type>::execute(sequence1, sequence2, edit_distance_cfg);
               };
    }

    template <typename is_semi_global_t>
    struct algorithm_impl
    {
        template <typename sequence1_t, typename sequence2_t, typename edit_distance_cfg_t>
        static auto
        execute(sequence1_t const & sequence1, sequence2_t const & sequence2, edit_distance_cfg_t edit_distance_cfg)
        {
            using alignment_result_t = alignment_result_type<sequence1_t, sequence2_t, edit_distance_cfg_t>;
            auto edit_distance_cfg_with_result_type =
                edit_distance_cfg | seqan3::align_cfg::detail::result_type<alignment_result_t>{};

            using edit_traits_t =
                seqan3::detail::default_edit_distance_trait_type<sequence1_t const &,
                                                                 sequence2_t const &,
                                                                 decltype(edit_distance_cfg_with_result_type),
                                                                 is_semi_global_t>;

            seqan3::detail::edit_distance_unbanded edit_distance{sequence1,
                                                                 sequence2,
                                                                 edit_distance_cfg_with_result_type,
                                                                 edit_traits_t{}};

            alignment_result_t align_result;
            edit_distance(0u,
                          [&align_result](auto && result)
                          {
                              align_result = std::forward<decltype(result)>(result);
                          });
            return align_result;
        }
    };
};

// ----------------------------------------------------------------------------
//  edit distance pairwise alignment (seqan2)
// ----------------------------------------------------------------------------

#ifdef SEQAN3_HAS_SEQAN2
struct edit_distance_algorithm_seqan2
{
    template <typename sequence1_t, typename sequence2_t, typename edit_distance_cfg_t>
    static auto select(edit_distance_cfg_t edit_distance_cfg)
    {
        using seqan3::get;

        // Note this converts in some cases a single configuration element like `seqan3::align_cfg::method_global`
        // to a seqan3::configuration to make seqan3::get<method_global> always work.
        seqan3::configuration align_cfg = edit_distance_cfg;
        auto method_global = get<seqan3::align_cfg::method_global>(align_cfg);
        assert(method_global.free_end_gaps_sequence2_leading == false);
        assert(method_global.free_end_gaps_sequence2_trailing == false);
        assert(method_global.free_end_gaps_sequence1_leading == method_global.free_end_gaps_sequence1_trailing);

        return method_global.free_end_gaps_sequence1_leading ?
               [] (sequence1_t & sequence1, sequence2_t & sequence2)
               {
                   return algorithm_impl<std::true_type>::execute(sequence1, sequence2);
               } :
               [] (sequence1_t & sequence1, sequence2_t & sequence2)
               {
                   return algorithm_impl<std::false_type>::execute(sequence1, sequence2);
               };
    }

    template <typename is_semi_global_t>
    struct algorithm_impl
    {
        template <typename sequence1_t, typename sequence2_t>
        static int execute(sequence1_t & sequence1, sequence2_t & sequence2)
        {
            constexpr bool is_semi_global = is_semi_global_t{};
            using method_t = std::conditional_t<is_semi_global, seqan2::MyersUkkonen, seqan2::MyersUkkonenGlobal>;
            using pattern_t = seqan2::Pattern<sequence2_t, method_t>;

            pattern_t pattern(sequence2, std::numeric_limits<int>::max());
            seqan2::Finder<sequence1_t> finder(sequence1);
            while (seqan2::find(finder, pattern))
                ;

            return -static_cast<int>(pattern.errors);
        }
    };
};
#endif // SEQAN3_HAS_SEQAN2

} // namespace seqan3::test
