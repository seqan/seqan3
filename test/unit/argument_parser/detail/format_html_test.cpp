// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <range/v3/view/remove_if.hpp>
#include <range/v3/algorithm/equal.hpp>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/argument_parser/detail/format_html.hpp>
#include <seqan3/io/stream/parse_condition.hpp>

using namespace seqan3;

TEST(html_test, html)
{
    std::string stdout;
    std::string expected;
    int option_value;
    bool flag_value;
    std::vector<std::string> pos_opt_value;

    // Empty html help page.
    const char * argv0[] = {"./help_add_test", "--export-help", "html"};
    argument_parser parser0("empty_options", 3, argv0);
    testing::internal::CaptureStdout();
    EXPECT_THROW(parser0.parse(), parser_interruption);
    stdout = testing::internal::GetCapturedStdout();
    expected = std::string("<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01//EN\" http://www.w3.org/TR/html4/strict.dtd\">"
                           "<html lang=\"en\">"
                           "<head>"
                           "<meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\">"
                           "<title> &mdash; </title>"
                           "</head>"
                           "<body>"
                           "<h1></h1>"
                           "<div></div>"
                           "<h2>Synopsis</h2>"
                           "<h2>Version</h2>"
                           "<strong>Last update:</strong> <br>"
                           "<strong> version:</strong> <br>"
                           "<strong>SeqAn version:</strong> 3.0.0<br>"
                           "</body></html>");
    EXPECT_TRUE(ranges::equal((stdout   | ranges::view::remove_if(is_space)),
                               expected | ranges::view::remove_if(is_space)));

   // Full html help page.
   argument_parser parser1("full", 3, argv0);
   parser1.info.synopsis.push_back("synopsis");
   parser1.info.synopsis.push_back("synopsis2");
   parser1.info.description.push_back("description");
   parser1.info.description.push_back("description2");
   parser1.info.short_description = "so short";
   parser1.info.url = "www.seqan.de";
   parser1.info.short_copyright = "short";
   parser1.info.long_copyright = "long";
   parser1.info.citation = "citation";
   parser1.add_option(option_value, 'i', "int", "this is a int option.");
   parser1.add_option(option_value, 'j', "jint", "this is a int option.");
   parser1.add_flag(flag_value, 'f', "flag", "this is a flag.");
   parser1.add_flag(flag_value, 'k', "kflag", "this is a flag.");
   parser1.add_positional_option(pos_opt_value, "this is a positional option.");
   parser1.add_positional_option(pos_opt_value, "this is a positional option.");
   parser1.info.examples.push_back("example");
   parser1.info.examples.push_back("example2");
   testing::internal::CaptureStdout();
   EXPECT_THROW(parser1.parse(), parser_interruption);
   stdout = testing::internal::GetCapturedStdout();
   expected = std::string("<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01//EN\" http://www.w3.org/TR/html4/strict.dtd\">"
                          "<html lang=\"en\">"
                          "<head>"
                          "<meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\">"
                          "<title> &mdash; </title>"
                          "</head>"
                          "<body>"
                          "<h1></h1>"
                          "<div></div>"
                          "<h2>Synopsis</h2>"
                          "<h2>Description</h2>"
                          "<p>"
                          "description"
                          "</p>"
                          "<p>"
                          "description2"
                          "</p>"
                          "<h2>Positional Arguments</h2>"
                          "<dl>"
                          "<dt><strong>ARGUMENT 2</strong> List of <em>STRING</em>'s</dt>"
                          "<dd>this is a positional option. </dd>"
                          "<dt><strong>ARGUMENT 2</strong> List of <em>STRING</em>'s</dt>"
                          "<dd>this is a positional option. </dd>"
                          "</dl>"
                          "<h2>Options</h2>"
                          "<dl>"
                          "<dt><strong>-i</strong>, <strong>--int</strong> <em>INT (32 bit)</em></dt>"
                          "<dd>this is a int option. </dd>"
                          "<dt><strong>-j</strong>, <strong>--jint</strong> <em>INT (32 bit)</em></dt>"
                          "<dd>this is a int option. </dd>"
                          "<dt><strong>-f</strong>, <strong>--flag</strong></dt>"
                          "<dd>this is a flag.</dd>"
                          "<dt><strong>-k</strong>, <strong>--kflag</strong></dt>"
                          "<dd>this is a flag.</dd>"
                          "</dl>"
                          "<h2>Version</h2>"
                          "<strong>Last update:</strong> <br>"
                          "<strong> version:</strong> <br>"
                          "<strong>SeqAn version:</strong> 3.0.0<br>"
                          "</body></html>");
   EXPECT_TRUE(ranges::equal((stdout   | ranges::view::remove_if(is_space)),
                              expected | ranges::view::remove_if(is_space)));
}
