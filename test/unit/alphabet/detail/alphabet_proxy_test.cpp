// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <seqan3/alphabet/detail/alphabet_proxy.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

#include "../alphabet_test_template.hpp"
#include "../alphabet_constexpr_test_template.hpp"

using namespace seqan3;

class alphabet_proxy_example : public alphabet_proxy<alphabet_proxy_example, dna4>
{
private:
    using base_t = alphabet_proxy<alphabet_proxy_example, dna4>;
    friend base_t;

    constexpr void on_update() noexcept
    {}

public:
    constexpr alphabet_proxy_example() noexcept : base_t{} {}
    constexpr alphabet_proxy_example(alphabet_proxy_example const &) = default;
    constexpr alphabet_proxy_example(alphabet_proxy_example &&) = default;
    constexpr alphabet_proxy_example & operator=(alphabet_proxy_example const &) = default;
    constexpr alphabet_proxy_example & operator=(alphabet_proxy_example &&) = default;
    ~alphabet_proxy_example() = default;
};

INSTANTIATE_TYPED_TEST_CASE_P(alphabet_proxy, alphabet, alphabet_proxy_example);
INSTANTIATE_TYPED_TEST_CASE_P(alphabet_proxy, alphabet_constexpr, alphabet_proxy_example);
