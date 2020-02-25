// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides cereal functionality for tests.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <gtest/gtest.h>

#include <fstream>

#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/core/platform.hpp>
#include <seqan3/core/type_traits/basic.hpp>
#include <seqan3/test/tmp_filename.hpp>

#if SEQAN3_WITH_CEREAL
#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/archives/xml.hpp>
#include <cereal/types/vector.hpp>
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
 * \param l The object to cerealise.
 */
template <cereal_input_archive in_archive_t, cereal_output_archive out_archive_t, typename value_t>
//!\cond
    requires cerealisable<value_t, in_archive_t, out_archive_t>
//!\endcond
void do_cerealisation(value_t && l)
{
    tmp_filename filename{"cereal_test"};

    {
        std::ofstream os{filename.get_path(), std::ios::binary};
        out_archive_t oarchive{os};
        oarchive(l);
    }

    {
        remove_cvref_t<value_t> in_l{};
        std::ifstream is{filename.get_path(), std::ios::binary};
        in_archive_t iarchive{is};
        iarchive(in_l);
        EXPECT_EQ(l, in_l);
    }
}

/*!\brief Tests if an object is serialise for all cereal archive types.
 * \tparam value_t The type to serialise.
 * \param l The object to serialise.
 *
 * If cereal is **not** available, this function is a NOP.
 * Otherwise it will call do_cerealisation() with cereal's `Binary`, `PortableBinary`, `JSON` and `XML` archives.
 */
template <typename value_t>
void do_serialisation([[maybe_unused]] value_t && l)
{
#if SEQAN3_WITH_CEREAL
    do_cerealisation<cereal::BinaryInputArchive,         cereal::BinaryOutputArchive>        (l);
    do_cerealisation<cereal::PortableBinaryInputArchive, cereal::PortableBinaryOutputArchive>(l);
    do_cerealisation<cereal::JSONInputArchive,           cereal::JSONOutputArchive>          (l);
    do_cerealisation<cereal::XMLInputArchive,            cereal::XMLOutputArchive>           (l);
#endif // SEQAN3_WITH_CEREAL
}
//!\endcond

} // namespace seqan3::test

} //namespace seqan3
