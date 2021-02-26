#include <fstream>
#include <gtest/gtest.h>

#include <seqan3/test/tmp_directory.hpp>

TEST(snippet_tmp_directory, tmp_directory_)
{
    // create a directory folder
    seqan3::test::tmp_directory tmp{};

    // Some function that should creates temporary files and removes them again
    {
        std::ofstream ofs{tmp.path() / "somefile.txt"};
        ofs << "Hello World!";
        ofs.close();

        std::filesystem::remove(tmp.path() / "somefile.txt");
    }

    // check that everything was cleanup properly
    EXPECT_TRUE(tmp.empty());
}

