// SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides EXPECT_THROW_MSG.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <gtest/gtest.h>

#include <seqan3/core/platform.hpp>

#ifdef EXPECT_THROW_MSG
#    warning "EXPECT_THROW_MSG is already defined."
#else
#    define EXPECT_THROW_MSG(statement, expected_exception, expected_message)                                          \
        try                                                                                                            \
        {                                                                                                              \
            statement;                                                                                                 \
            std::string const message = "Expected: " #statement " throws an exception of type " #expected_exception    \
                                        ".\n  Actual: it throws nothing.";                                             \
            GTEST_NONFATAL_FAILURE_(message.data());                                                                   \
        }                                                                                                              \
        catch (expected_exception const & exception)                                                                   \
        {                                                                                                              \
            if (auto result = ::testing::internal::EqHelper::Compare("Expected",                                       \
                                                                     "Actual",                                         \
                                                                     std::string_view{expected_message},               \
                                                                     std::string_view{exception.what()});              \
                !result)                                                                                               \
            {                                                                                                          \
                std::string message = #statement " throws the correct exception, but the description is incorrect.\n"; \
                message += result.failure_message();                                                                   \
                GTEST_NONFATAL_FAILURE_(message.data());                                                               \
            }                                                                                                          \
        }                                                                                                              \
        catch (std::exception const & exception)                                                                       \
        {                                                                                                              \
            std::string message = "Expected: " #statement " throws an exception of type " #expected_exception ".\n  "; \
            message += "Actual: it throws ";                                                                           \
            message += ::testing::internal::GetTypeName(typeid(exception));                                            \
            message += " with description \"";                                                                         \
            message += exception.what();                                                                               \
            message += "\".";                                                                                          \
            GTEST_NONFATAL_FAILURE_(message.data());                                                                   \
        }                                                                                                              \
        catch (...)                                                                                                    \
        {                                                                                                              \
            std::string message = "Expected: " #statement " throws an exception of type " #expected_exception ".\n  "; \
            message += "Actual: it throws an unknown exception.";                                                      \
            GTEST_NONFATAL_FAILURE_(message.data());                                                                   \
        }
#endif
