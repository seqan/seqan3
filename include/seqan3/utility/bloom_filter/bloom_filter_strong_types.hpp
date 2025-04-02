// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides strong types for the (Interleaved) Bloom Filter.
 */

#pragma once

#include <seqan3/core/detail/strong_type.hpp>

//Todo: When removing search/dream_index/interleaved_bloom_filter.hpp, the contents of this header can be moved
// into utility/bloom_filter/bloom_filter.hpp

namespace seqan3
{

//!\brief Determines if the Interleaved Bloom Filter is compressed.
//!\ingroup utility_bloom_filter
enum data_layout : bool
{
    uncompressed, //!< The Interleaved Bloom Filter is uncompressed.
    compressed    //!< The Interleaved Bloom Filter is compressed.
};

//!\brief A strong type that represents the number of bins for the seqan3::interleaved_bloom_filter.
//!\ingroup utility_bloom_filter
struct bin_count : public detail::strong_type<size_t, bin_count, detail::strong_type_skill::convert>
{
    using detail::strong_type<size_t, bin_count, detail::strong_type_skill::convert>::strong_type;
};

//!\brief A strong type that represents the number of bits for each bin in the seqan3::interleaved_bloom_filter.
//!\ingroup utility_bloom_filter
struct bin_size : public detail::strong_type<size_t, bin_size, detail::strong_type_skill::convert>
{
    using detail::strong_type<size_t, bin_size, detail::strong_type_skill::convert>::strong_type;
};

//!\brief A strong type that represents the number of hash functions for the seqan3::interleaved_bloom_filter.
//!\ingroup utility_bloom_filter
struct hash_function_count : public detail::strong_type<size_t, hash_function_count, detail::strong_type_skill::convert>
{
    using detail::strong_type<size_t, hash_function_count, detail::strong_type_skill::convert>::strong_type;
};

//!\brief A strong type that represents the bin index for the seqan3::interleaved_bloom_filter.
//!\ingroup utility_bloom_filter
struct bin_index : public detail::strong_type<size_t, bin_index, detail::strong_type_skill::convert>
{
    using detail::strong_type<size_t, bin_index, detail::strong_type_skill::convert>::strong_type;
};

} // namespace seqan3
