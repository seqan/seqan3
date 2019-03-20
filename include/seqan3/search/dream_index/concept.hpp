// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides concepts for the seqan3::dream_index.
 */

#pragma once

#include <seqan3/core/metafunction/template_inspection.hpp>
#include <seqan3/range/container/concept.hpp>
#include <seqan3/search/dream_index/binning_directory.hpp>
#include <seqan3/search/dream_index/detail/bitvector.hpp>

namespace seqan3
{

    /*!\interface seqan3::DreamIndexTraits <>
     * \brief Concept for DREAM index traits.
     *
     * This concept defines the traits for the DREAM index.
     */
    //!\cond
    template <typename t>
    SEQAN3_CONCEPT DreamIndexTraits = requires (t v)
    {
        typename t::alphabet_t;
        typename t::bitvector_strategy;
        typename t::directory_strategy;

        requires Alphabet<typename t::alphabet_t>;
        requires reservable_container_concept<detail::bitvector<typename t::bitvector_strategy>>;
        // requires detail::is_type_specialisation_of_v<binning_directory<typename t::directory_strategy>,
        //                                              binning_directory>;
        requires std::Same<typename t::directory_strategy, direct> || std::Same<typename t::directory_strategy, ibf>;
    };
    //!\endcond
    /*!\name Requirements for seqan3::DreamIndexTraits
     * \relates seqan3::DreamIndexTraits
     * \brief The DREAM index traits must provide the following types:
     * \{
     *
     * \typedef typename typename t::alphabet_t alphabet_t
     * \brief Declares the alphabet type for building the index.
     *        Must satisfy seqan3::Alphabet.
     *
     * \typedef typename typename t::bitvector_strategy bitvector_strategy
     * \brief Declares the type of the underlying bitvector structure.
     *        The resulting seqan3::detail::bitvector specialisation must satisfy seqan3::reservable_container_concept.
     *
     * \typedef typename typename t::directory_strategy directory_strategy
     * \brief Declares the type of the underlying binning directory.
     *        The resulting seqan3::binning_directory specialisation must satisfy seqan3::BinningDirectory.
     *
     * \}
     */

} // namespace seqan3
