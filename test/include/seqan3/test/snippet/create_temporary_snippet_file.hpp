// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <seqan3/std/filesystem>
#include <fstream>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/test/tmp_directory.hpp>

//!\cond
struct create_temporary_snippet_file
{
    std::filesystem::path file_path;

    create_temporary_snippet_file(std::filesystem::path const & file_name, std::string const & file_raw)
        : file_path{}
    {
        // create single folder instance (across multiple create_temporary_snippet_file instances) that is valid during
        // the complete program. std::filesystem::current_path will point to that location.
        static seqan3::test::tmp_directory const tmp_folder{[]()
        {
            seqan3::test::tmp_directory tmp{};
            std::filesystem::current_path(tmp.path());
            return tmp;
        }()};

        file_path = tmp_folder.path();
        file_path /= file_name;

        // create file if exists
        if (!file_raw.empty())
        {
            std::ofstream file{file_path};
            file << file_raw.substr(1); // skip first newline
        }
    }

    ~create_temporary_snippet_file()
    {
        std::error_code ec{};
        std::filesystem::remove(file_path, ec);

        if (ec)
            seqan3::debug_stream << "[WARNING] Could not delete " << file_path << ". " << ec.message() << '\n';
    }
};
//!\endcond
