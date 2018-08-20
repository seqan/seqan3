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
// ============================================================================

#include <gtest/gtest.h>

#include <seqan3/alignment/configuration/utility.hpp>
#include <seqan3/core/algorithm/config_element_access.hpp>
#include <seqan3/core/algorithm/config_element_base.hpp>
#include <seqan3/core/algorithm/configuration.hpp>

using namespace seqan3;

TEST(utility, align_cfg_id)
{
    // NOTE(rrahn): You must update this test if you add a new value to align_cfg::id
    EXPECT_EQ(static_cast<uint8_t>(align_cfg::id::SIZE), 1);
}

class bar : public detail::config_element_base<bar>
{
    friend class detail::config_element_access<bar>;

    int value;
};

TEST(utility, on_align_config)
{
    EXPECT_TRUE((std::is_same_v<typename detail::on_align_config<align_cfg::id::SIZE>::invoke<bar>, std::false_type>));
}

TEST(utility, align_config_type_to_id)
{
    EXPECT_EQ(detail::align_config_type_to_id<bar>::value, align_cfg::id::SIZE);
    EXPECT_EQ(detail::align_config_type_to_id_v<bar>, align_cfg::id::SIZE);
}
