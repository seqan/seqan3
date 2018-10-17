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

/*!\file
 * \brief This header includes C++17 filesystem support and imports it into namespace seqan3::filesystem (independent of whether it is marked as "experimental").
 * \author Vitor C. Piro <pirov AT zedat.fu-berlin.de >
 */

#pragma once


// Temporal workaround GCC 8.1 not implementing filesystem for windows
#ifndef BOOST_FILESYSTEM_FORCE


#if __has_include(<filesystem>)
#include <filesystem>
#else
#include <experimental/filesystem>
#endif // __has_include(experimental/filesystem)

#include <seqan3/core/platform.hpp>

namespace seqan3
{
#if __has_include(<filesystem>)
namespace filesystem = std::filesystem;
#else
namespace filesystem = std::experimental::filesystem;
#endif // __has_include(experimental/filesystem)
} // namespace seqan3


#else    // defined BOOST_FILESYSTEM_FORCE
#define _GLIBCXX_FILESYSTEM 1
#  include <chrono>

// Workaround boost assert producing warnings (at least in 1.67):
//   C:/MinGW/include/boost/mpl/assert.hpp:188:21: error: unnecessary parentheses in declaration of 'assert_arg' [-Werror=parentheses]
//    failed ************ (Pred::************
//                     ^
//   C:/MinGW/include/boost/mpl/assert.hpp:193:21: error: unnecessary parentheses in declaration of 'assert_not_arg' [-Werror=parentheses]
//    failed ************ (boost::mpl::not_<Pred>::************
// todo: report to boost?? suppress error only for those two lines
#pragma GCC diagnostic push
#  pragma GCC diagnostic warning "-Wparentheses"
#  include <boost/filesystem.hpp>
#pragma GCC diagnostic pop
namespace std {
    namespace filesystem {
        using namespace boost::filesystem;
        using file_time_type = std::chrono::time_point<std::chrono::system_clock>;

        enum class file_type {
            none      = boost::filesystem::file_type::status_unknown,
            not_found = boost::filesystem::file_type::file_not_found,
            regular   = boost::filesystem::file_type::regular_file,
            directory = boost::filesystem::file_type::directory_file,
            symlink   = boost::filesystem::file_type::symlink_file,
            block     = boost::filesystem::file_type::block_file,
            character = boost::filesystem::file_type::character_file,
            fifo      = boost::filesystem::file_type::fifo_file,
            socket    = boost::filesystem::file_type::socket_file,
            unknown   = boost::filesystem::file_type::type_unknown
        };
    } // filesystem
      // std::string operator ()(const filesystem::path& p){return p.string();}
} // std
namespace seqan3 {
    namespace filesystem = std::filesystem;
}

#endif //   BOOST_FILESYSTEM_FORCE




