// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>
#include <fstream>

#include <seqan3/test/tmp_filename.hpp>

namespace fs = std::filesystem;

// aggregate initialization
TEST(tmp_filename_aggr, aggr)
{
    seqan3::test::tmp_filename t1{"aggr_test"};
    seqan3::test::tmp_filename t2("aggr_test");
    EXPECT_NE(t1.get_path(), t2.get_path());
    EXPECT_TRUE(fs::exists(t1.get_path().parent_path()));
    EXPECT_TRUE(fs::exists(t2.get_path().parent_path()));

    EXPECT_TRUE(fs::equivalent(fs::temp_directory_path(), t1.get_path().parent_path().parent_path()));
    EXPECT_TRUE(fs::equivalent(fs::temp_directory_path(), t2.get_path().parent_path().parent_path()));
}

// nullptr as filename
TEST(tmp_filename_nullptr, null_ptr)
{
    EXPECT_THROW(seqan3::test::tmp_filename t1{nullptr}, fs::filesystem_error);
    EXPECT_THROW(seqan3::test::tmp_filename t1(nullptr), fs::filesystem_error);
}

// move construction
TEST(tmp_filename_mv_ctr, mv_ctr)
{
    seqan3::test::tmp_filename t1{"mv_ctr_test"};
    seqan3::test::tmp_filename t2{"mv_ctr_test"};
    seqan3::test::tmp_filename t3{std::move(t2)};
    EXPECT_NE(t1.get_path(), t3.get_path());
    seqan3::test::tmp_filename t4(std::move(t1));
    EXPECT_NE(t3.get_path(), t4.get_path());
}

// move assignment
TEST(tmp_filename_mv_assign, mv_assign)
{
    seqan3::test::tmp_filename t1{"mv_ctr_test"};
    seqan3::test::tmp_filename t2{"mv_ctr_test"};
    seqan3::test::tmp_filename t3 = std::move(t2);
    EXPECT_NE(t1.get_path(), t3.get_path());
}

// destructor
TEST(tmp_filename_dtr, dtr)
{
    auto t1 = std::make_unique<seqan3::test::tmp_filename>("delete_test");
    auto path = t1->get_path();
    std::ofstream os{path, std::ios::out};
    os << "delete_test";
    os.close();
    EXPECT_TRUE(fs::exists(path));
    EXPECT_TRUE(fs::exists(path.parent_path()));
    t1.reset();
    EXPECT_FALSE(fs::exists(path));
    EXPECT_FALSE(fs::exists(path.parent_path()));
}

//!\brief A wrapper class to handle construction, destruction and permission handling of a temporary directory.
class read_only_directory
{
public:
    read_only_directory() = delete; //< Deleted.
    read_only_directory(read_only_directory const &) = delete; //< Deleted.
    read_only_directory(read_only_directory &&) = default; //< Defaulted.
    read_only_directory & operator=(read_only_directory const &) = delete; //< Deleted.
    read_only_directory & operator=(read_only_directory &&) = default; //< Defaulted.

    /*!\brief Construct a directory from a given path.
     * \throws if std::filesystem::create_directory throws.
     * \throws std::filesystem::filesystem_error if the directory already exists.
     * \throws if std::filesystem::permissions throws.
     */
    explicit read_only_directory(std::filesystem::path path) : directory_path(std::move(path))
    {
        if (!std::filesystem::create_directory(directory_path))
        {
            throw std::filesystem::filesystem_error{"The read_only_directory path already exists",
                                                    directory_path,
                                                    std::make_error_code(std::errc::file_exists)};
        }
        std::filesystem::permissions(directory_path,
                                     std::filesystem::perms::owner_write,
                                     std::filesystem::perm_options::remove);
    }

    //!\brief Destructor tries to delete the directory.
    ~read_only_directory()
    {
        [[maybe_unused]] std::error_code ec;
        release_impl(ec);
    }

    /*!\brief Tries to remove the directory.
     * \throws if std::filesystem::permissions throws.
     * \throws if std::filesystem::remove throws.
     * \throws std::filesystem::filesystem_error if the directory does not exist.
     */
    void release() const
    {
        std::error_code ec;
        release_impl(ec);

        if (ec.value())
            throw std::filesystem::filesystem_error("Error while removing directory", directory_path, ec);
    }

private:
    //!/brief The path of the managed directory.
    std::filesystem::path const directory_path;

    /*!\brief Tries to remove the directory.
     */
    void release_impl(std::error_code & ec) const noexcept
    {
        std::filesystem::permissions(directory_path,
                                     std::filesystem::perms::owner_write,
                                     std::filesystem::perm_options::add,
                                     ec);
        std::filesystem::remove(directory_path, ec);
    }
};

TEST(tmp_filename_throw, directory_not_writeable)
{
    try
    {
        // Create a directory in the temporary directory.
        std::filesystem::path test_path = std::filesystem::temp_directory_path();
        test_path /= "seqan3_tmp_filename_throw";
        read_only_directory test_directory{test_path};

        // Set TMPDIR. This is the first env var that is looked up for `temp_directory_path` inside the `tmp_filename`.
        setenv("TMPDIR", test_path.c_str(), 1); // name, value, overwrite

        EXPECT_THROW(seqan3::test::tmp_filename t1{"throw"}, std::filesystem::filesystem_error);

        // Remove directory.
        test_directory.release();
    }
    catch (std::exception const & e)
    {
        FAIL() << e.what();
    }
}
