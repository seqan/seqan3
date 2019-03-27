// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Tobias Loka <LokaT AT rki.de>
 * \brief Provides functions to test the fulfillment of range concepts.
 */

#include <gtest/gtest.h>

#include <seqan3/io/stream/debug_stream.hpp>
#include <seqan3/std/ranges>

using namespace seqan3;

namespace seqan3::test 
{

//!\brief Enum for the different types of range concepts.
enum ConceptType : uint8_t
{
    Input,
    Forward,
    Bidirectional,
    RandomAccess,
    Contiguous,
    Viewable,
    View,
    Sized,
    Common,
    Output,
    ConstIterable
};

/*!\brief Convert a ConceptType object to std::string for output.
 * \tparam conc The range concept to convert to string as seqan3::test::ConceptType object.
 * \return The concept type as std::string.
 *
 *\details
 * 
 * The return string is not a 1:1 conversion of the enum type name but rather the full name
 * of the range concept it stands for.
 */
std::string to_string(ConceptType const conc)
{
    switch(conc)
    {
        case Input:
            return "InputRange";
        case Forward:
            return "ForwardRange";
        case Bidirectional:
            return "BidirectionalRange";
        case RandomAccess:
            return "RandomAccessRange";
        case Contiguous:
            return "ContiguousRange";
        case Viewable:
            return "ViewableRange";
        case View:
            return "View";
        case Sized:
            return "SizedRange";
        case Common:
            return "CommonRange";
        case Output:
            return "OutputRange";
        case ConstIterable:
            return "const_iterable_concept";
        default:
            return "Undefined";
    }
}

/*!\brief Expression to check the fulfillment of a range concept.
 * \tparam conc    The range concept to check as seqan3::test::ConceptType object.
 * \return True if TYPE fulfills conc. False if not.
 * \throws std::logic_error    Throws logic error if the concept type given (concept) is unknown.
 *
 * ### Example
 * using namespace seqan3::test;
 * std::string vec{"HelloWorld"};
 * EXPECT_TRUE(fulfilled<decltype(vec)>(ConceptType::Input));
 */
template<typename TYPE>
constexpr bool fulfilled(ConceptType const conc)
{
    switch (conc)
    {
        case Input:
            return std::ranges::InputRange<TYPE>;
        case Forward:
            return std::ranges::ForwardRange<TYPE>;
        case Bidirectional:
            return std::ranges::BidirectionalRange<TYPE>;
        case RandomAccess:
            return std::ranges::RandomAccessRange<TYPE>;
        case Contiguous:
            return std::ranges::ContiguousRange<TYPE>;
        case Viewable:
            return std::ranges::ViewableRange<TYPE>;
        case View:
            return std::ranges::View<TYPE>;
        case Sized:
            return std::ranges::SizedRange<TYPE>;
        case Common:
            return std::ranges::CommonRange<TYPE>;
        case Output:
            return std::ranges::OutputRange<TYPE, typename seqan3::value_type<TYPE>::type>;
        case ConstIterable:
            return const_iterable_concept<TYPE>;
        default:
            throw std::logic_error("Unexpected ConceptType in view_concept_check::fulfilled(...).");
    }
}

/*!\brief Checks if the fulfillment of a concept is preserved from IN_TYPE to OUT_TYPE.
 * \tparam concepts    A list of range concepts to check as seqan3::test::ConceptType object.
 * \return True if fulfillment of conc is preserved from IN_TYPE to OUT_TYPE. False if not.
 *
 * ### Example
 * using namespace seqan3::test;
 * std::string vec{"HelloWorld"};
 * auto view = vec | std::view::reverse;
 * EXPECT_TRUE((preserved<decltype(vec),decltype(view)>({Input, Forward}))); //double brackets required
 *
 *\details
 *
 * This function is designed for use inside EXPECT_TRUE(...). Using other test scenarios (e.g., EXPECT_FALSE)
 * can lead to unintended command line output via debug_stream.
 */
template<std::ranges::Range IN_TYPE, typename OUT_TYPE>
constexpr bool preserved(std::initializer_list<ConceptType> const & concepts)
{
    bool success = true;
    auto it = concepts.begin();
    for (; it != concepts.end(); ++it)
    {
        if(fulfilled<IN_TYPE>(*it) != fulfilled<OUT_TYPE>(*it))
        {
            debug_stream << "Preserved check of concept '" << to_string(*it) << "' failed:\n" 
                         << " IN_TYPE: " << fulfilled<IN_TYPE>(*it) << '\n'
                         << "OUT_TYPE: " << fulfilled<OUT_TYPE>(*it) << '\n'
                         << "are expected to be equal.\n";
            success = false;
        }
    }
    return success;
}

/*!\brief Checks if the fulfillment of a concept is lost from IN_TYPE to OUT_TYPE.
 * \tparam concepts    A list of range concepts to check as seqan3::test::ConceptType object.
 * \return True if fulfillment of conc is lost from IN_TYPE to OUT_TYPE. False if not.
 *
 * ### Example
 * using namespace seqan3::test;
 * std::string vec{"HelloWorld"};
 * auto view = vec | std::view::reverse;
 * EXPECT_TRUE((lost<decltype(vec), decltype(view)>({Contiguous}))); // double brackets required
 *
 * \details
 *
 * The lost function checks a real loss, i.e. true for IN_TYPE and false for OUT_TYPE. 
 * You can use weak_lost(...) function to only check false for OUT_TYPE.
 * This function is designed for use inside EXPECT_TRUE(...). Using other test scenarios (e.g., EXPECT_FALSE)
 * can lead to unintended command line output via debug_stream.
 */
template<std::ranges::Range IN_TYPE, typename OUT_TYPE>
constexpr bool lost(std::initializer_list<ConceptType> const & concepts)
{
    bool success = true;
    auto it = concepts.begin();
    for (; it != concepts.end(); ++it)
    {
        if(!fulfilled<IN_TYPE>(*it) || fulfilled<OUT_TYPE>(*it))
        {
            debug_stream << "Lost check of concept '" << to_string(*it) << "' failed:\n" 
                         << " IN_TYPE: " << fulfilled<IN_TYPE>(*it) << " (expected true)" << '\n'
                         << "OUT_TYPE: " << fulfilled<OUT_TYPE>(*it) << " (expected false)" <<'\n';
            success = false;
        }
    }
    return success;
}

/*!\brief Checks weakly if the fulfillment of a concept is lost from IN_TYPE to OUT_TYPE.
 * \tparam concepts    A list of range concepts to check as seqan3::test::ConceptType object.
 * \return True if fulfillment of conc is (weakly) lost in OUT_TYPE. False if not.

 *
 * ### Example
 * using namespace seqan3::test;
 * std::string vec{"HelloWorld"};
 * auto view = vec | std::view::reverse;
 * EXPECT_TRUE(weak_lost<decltype(view)>({Contiguous}));
 *
 * \details
 *
 * The lost function weakly checks a loss, i.e. only false for OUT__TYPE. 
 * You can use lost(...) function to check true for IN_TYPE and false for OUT_TYPE.
 * This function is designed for use inside EXPECT_TRUE(...). Using other test scenarios (e.g., EXPECT_FALSE)
 * can lead to unintended command line output via debug_stream.
 */
template<typename OUT_TYPE>
constexpr bool weak_lost(std::initializer_list<ConceptType> const & concepts)
{
    bool success = true;
    auto it = concepts.begin();
    for (; it != concepts.end(); ++it)
    {
        if(fulfilled<OUT_TYPE>(*it))
        {
            debug_stream << "Weak lost check of concept '" << to_string(*it) << "' failed:\n" 
                         << "OUT_TYPE: " << fulfilled<OUT_TYPE>(*it) << " (expected false)" <<'\n';
            success = false;
        }
    }
    return success;
}

/*!\brief Checks if the fulfillment of a concept is guaranteed from IN_TYPE to OUT_TYPE.
 * \tparam concepts    A list of range concepts to check as seqan3::test::ConceptType object.
 * \return True if fulfillment of conc is guaranteed in OUT_TYPE. False if not.
 *
 * ### Example
 * using namespace seqan3::test;
 * std::string vec{"HelloWorld"};
 * auto view = vec | std::view::reverse;
 * // This EXPECT_TRUE will fail, as InputRange concept is fulfilled by IN_TYPE.
 * // For this setup, you would have to use weak_guaranteed(...) function instead.
 * EXPECT_TRUE((guaranteed<decltype(vec),decltype(view)>({Input}))); // double brackets required
 *
 * \details
 *
 * The lost function checks a real guarantee, i.e. false for IN_TYPE and true for OUT_TYPE. 
 * You can use weak_guarantee(...) function to only check true for OUT_TYPE.
 * This function is designed for use inside EXPECT_TRUE(...). Using other test scenarios (e.g., EXPECT_FALSE)
 * can lead to unintended command line output via debug_stream.
 */
template<typename IN_TYPE, typename OUT_TYPE>
constexpr bool guaranteed(std::initializer_list<ConceptType> const & concepts)
{
    bool success = true;
    auto it = concepts.begin();
    for (; it != concepts.end(); ++it)
    {
        if(fulfilled<IN_TYPE>(*it) || !fulfilled<OUT_TYPE>(*it))
        {
            debug_stream << "Guaranteed check of concept '" << to_string(*it) << "' failed:\n" 
                         << " IN_TYPE: " << fulfilled<IN_TYPE>(*it) << " (expected false)" << '\n'
                         << "OUT_TYPE: " << fulfilled<OUT_TYPE>(*it) << " (expected true)" <<'\n';
            success = false;
        }
    }
    return success;
}

/*!\brief Checks weakly if the fulfillment of a concept is lost from IN_TYPE to OUT_TYPE.
 * \tparam concepts    A list of range concepts to check as seqan3::test::ConceptType object.
 * \return True if fulfillment of conc is (weakly) guaranteed in OUT_TYPE. False if not.
 *
 * ### Example
 * using namespace seqan3::test;
 * std::string vec{"HelloWorld"};
 * auto view = vec | std::view::reverse;
 * EXPECT_TRUE(weak_guaranteed<decltype(view)>({Input}));
 *
 * \details
 *
 * The weak_guaranteed(...) function weakly checks a guarantee, i.e. only true for OUT_TYPE. 
 * You can use guaranteed(...) function to check true for IN_TYPE and true for OUT_TYPE.
 * This function is designed for use inside EXPECT_TRUE(...). Using other test scenarios (e.g., EXPECT_FALSE)
 * can lead to unintended command line output via debug_stream.
 */
template<typename OUT_TYPE>
constexpr bool weak_guaranteed(std::initializer_list<ConceptType> const & concepts)
{
    bool success = true;
    auto it = concepts.begin();
    for (; it != concepts.end(); ++it)
    {
        if(!fulfilled<OUT_TYPE>(*it))
        {
            debug_stream << "Weak guaranteed check of concept '" << to_string(*it) << "' failed:\n" 
                         << "OUT_TYPE: " << fulfilled<OUT_TYPE>(*it) << " (expected true)" <<'\n';
            success = false;
        }
    }
    return success;
}
}