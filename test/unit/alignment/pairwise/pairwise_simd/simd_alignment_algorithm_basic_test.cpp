// // -----------------------------------------------------------------------------------------------------
// // Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// // Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// // This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// // shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// // -----------------------------------------------------------------------------------------------------
//
// #include <gtest/gtest.h>
//
// #include <seqan3/alphabet/nucleotide/dna4.hpp>
// #include <seqan3/alignment/configuration/all.hpp>
// #include <seqan3/alignment/pairwise/policy/all.hpp>
// #include <seqan3/alignment/scoring/simd_gap_scheme.hpp>
// #include <seqan3/alignment/scoring/simd_scoring_scheme_simple.hpp>
// #include <seqan3/range/container/aligned_allocator.hpp>
//
// using namespace seqan3;
//
// TEST(basic_simd_alignment_test, global_unbanded_score_only_equal_length)
// {
//
//     using simd_t = simd::simd_type_t<int32_t>;
//
//     using cell_t = std::tuple<simd_t, simd_t, detail::ignore_t>;
//     using allocator_t = aligned_allocator<cell_t, sizeof(simd_t)>;
//     using affine_gap_policy_t = detail::deferred_crtp_base<detail::affine_gap_policy, cell_t>;
//     using affine_gap_init_policy_t = detail::deferred_crtp_base<detail::affine_gap_init_policy>;
//     using matrix_policy_t = detail::deferred_crtp_base<detail::unbanded_score_dp_matrix_policy, allocator_t>;
//     using find_optimum_t = detail::deferred_crtp_base<detail::find_optimum_policy>;
//
//     // Now we need an align_config that can deal with the score.
//
//     auto config = align_cfg::mode{global_alignment} |
//                   align_cfg::scoring{simd_scoring_scheme_simple<simd_t>{match_score{4}, mismatch_score{-5}}} |
//                   align_cfg::gap{simd_gap_scheme<simd_t>{gap_score{-1}, gap_open_score{-11}}};
//
//     using alignment_algorithm_t = detail::alignment_algorithm<decltype(config), affine_gap_policy_t, affine_gap_init_policy_t, matrix_gap_policy_t, find_optimum_t>;
//
//     [[maybe_unused]] alignment_algorithm_t algo{config};
//
//     // auto seq1 = "ACGTACGATCGAC"_dna4;
//     // auto seq2 = "AGATCGAGCTACGAGCTA"_dna4;
//     //
//     // using native_simd_t = simd_type_t<int16_t>;
//     // constexpr size_t v = simd_traits<native_simd_t>::length;
//     //
//     // // how can we make this programatically.
//     //
//     // iter fwd_iter = begin(fwd_range<sequences>)
//     // iter_vec;
//     // while (in_avail() && fwd_iter != end(fwd_iter<sequences>))
//     // {
//     //     iter_vec.push_back(std::pair{begin(*fwd_iter), end(*fwd_iter)});
//     // }
//     //
//     // range {range} // max > simd_size
//     //
//     // auto soa_vec = transform_to_soa(range_ranges, padding);
//
//
//
//
//
// }
