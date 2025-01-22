// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alignment/pairwise/alignment_configurator.hpp>

int main()
{
    using sequences_t = std::vector<std::pair<std::string, std::string>>;
    using config_t = decltype(seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme
                              | seqan3::align_cfg::output_score{});

    using first_seq_t = std::tuple_element_t<0, std::ranges::range_value_t<sequences_t>>;
    using second_seq_t = std::tuple_element_t<1, std::ranges::range_value_t<sequences_t>>;

    // Select the result type based on the sequences and the configuration.
    using result_t =
        seqan3::alignment_result<typename seqan3::detail::align_result_selector<std::remove_reference_t<first_seq_t>,
                                                                                std::remove_reference_t<second_seq_t>,
                                                                                config_t>::type>;
    // Define the function wrapper type.
    using function_wrapper_t = std::function<result_t(first_seq_t &, second_seq_t &)>;

    static_assert(seqan3::detail::is_type_specialisation_of_v<function_wrapper_t, std::function>);
}
