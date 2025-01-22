// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides cereal functionality for tests.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <gtest/gtest.h>

#include <fstream>

#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/core/platform.hpp>
#include <seqan3/test/tmp_directory.hpp>
#include <seqan3/utility/type_traits/basic.hpp>

#if SEQAN3_WITH_CEREAL
#    include <cereal/archives/binary.hpp>
#    include <cereal/archives/json.hpp>
#    include <cereal/archives/portable_binary.hpp>
#    include <cereal/archives/xml.hpp>
#    include <cereal/types/vector.hpp>
#endif // SEQAN3_WITH_CEREAL

namespace seqan3
{

namespace test
{
//!\cond DEV
/*!\brief Tests if an object is cerealisable.
 * \tparam in_archive_t  Type of the cereal input archive. Must model seqan3::cereal_input_archive.
 * \tparam out_archive_t Type of the cereal output archive. Must model seqan3::cereal_output_archive.
 * \tparam value_t       The type to cerealise. Must model seqan3::cerealisable.
 * \param value The object to cerealise.
 */
template <cereal_input_archive in_archive_t, cereal_output_archive out_archive_t, typename value_t>
    requires cerealisable<value_t, in_archive_t, out_archive_t>
void do_cerealisation(value_t && value)
{
    tmp_directory tmp{};
    auto filename = tmp.path() / "cereal_test";

    {
        std::ofstream os{filename, std::ios::binary};
        out_archive_t oarchive{os};
        oarchive(value);
    }

    {
        std::remove_cvref_t<value_t> value_from_archive{};
        std::ifstream is{filename, std::ios::binary};
        in_archive_t iarchive{is};
        iarchive(value_from_archive);
        EXPECT_TRUE(value == value_from_archive);
    }
}

/*!\brief Tests if an object is serialise for all cereal archive types.
 * \tparam value_t The type to serialise.
 * \param value The object to serialise.
 *
 * If cereal is **not** available, this function is a NOP.
 * Otherwise it will call do_cerealisation() with cereal's `Binary`, `PortableBinary`, `JSON` and `XML` archives.
 */
template <typename value_t>
void do_serialisation([[maybe_unused]] value_t && value)
{
#if SEQAN3_WITH_CEREAL
    do_cerealisation<cereal::BinaryInputArchive, cereal::BinaryOutputArchive>(value);
    do_cerealisation<cereal::PortableBinaryInputArchive, cereal::PortableBinaryOutputArchive>(value);
    do_cerealisation<cereal::JSONInputArchive, cereal::JSONOutputArchive>(value);
    do_cerealisation<cereal::XMLInputArchive, cereal::XMLOutputArchive>(value);
#endif // SEQAN3_WITH_CEREAL
}
//!\endcond

} // namespace test

} //namespace seqan3
