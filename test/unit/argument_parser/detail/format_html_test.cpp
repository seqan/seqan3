// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <range/v3/view/remove_if.hpp>
#include <range/v3/algorithm/equal.hpp>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/argument_parser/detail/format_html.hpp>
#include <seqan3/io/stream/parse_condition.hpp>

using namespace seqan3;

TEST(html_format, empty_information)
{
    std::string my_stdout;
    std::string expected;

    // Empty html help page.
    const char * argv0[] = {"./help_add_test", "--export-help", "html"};
    argument_parser parser0{"empty_options", 3, argv0};
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser0.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    my_stdout = testing::internal::GetCapturedStdout();
    expected = std::string("<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01//EN\" http://www.w3.org/TR/html4/strict.dtd\">\n"
                           "<html lang=\"en\">\n"
                           "<head>\n"
                           "<meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\">\n"
                           "<title>empty_options &mdash; </title>\n"
                           "</head>\n"
                           "<body>\n"
                           "<h1>empty_options</h1>\n"
                           "<div></div>\n"
                           "<h2>Options</h2>\n"
                           "<h3>Basic options:</h3>\n"
                           "<dl>\n"
                           "<dt><strong>-h</strong>, <strong>--help</strong></dt>\n"
                           "<dd>Prints the help page.</dd>\n"
                           "<dt><strong>-hh</strong>, <strong>--advanced-help</strong></dt>\n"
                           "<dd>Prints the help page including advanced options.</dd>\n"
                           "<dt><strong>--version</strong></dt>\n"
                           "<dd>Prints the version information.</dd>\n"
                           "<dt><strong>--copyright</strong></dt>\n"
                           "<dd>Prints the copyright/license information.</dd>\n"
                           "<dt><strong>--export-help</strong> (std::string)</dt>\n"
                           "<dd>Export the help page information. Value must be one of [html, man].</dd>\n"
                           "</dl>\n"
                           "<h3></h3>\n"
                           "<h2>Version</h2>\n"
                           "<strong>Last update:</strong> <br>\n"
                           "<strong>empty_options version:</strong> <br>\n"
                           "<strong>SeqAn version:</strong> 3.0.0<br>\n"
                           "<br>\n"
                           "</body></html>");
    EXPECT_EQ(my_stdout, expected);

    const char * argv1[] = {"./help_add_test", "--export-help=html"};
    argument_parser parser1{"empty_options", 2, argv1};
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser1.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    my_stdout = testing::internal::GetCapturedStdout();
    EXPECT_EQ(my_stdout, expected);
}

TEST(html_format, full_information_information)
{
    std::string my_stdout;
    std::string expected;
    int option_value{5};
    bool flag_value;
    std::vector<std::string> pos_opt_value{};

   // Full html help page.
   const char * argv0[] = {"./help_add_test", "--export-help", "html"};
   argument_parser parser1{"program_full_options", 3, argv0};
   parser1.info.synopsis.push_back("./some_binary_name synopsis");
   parser1.info.synopsis.push_back("./some_binary_name synopsis2");
   parser1.info.description.push_back("description");
   parser1.info.description.push_back("description2");
   parser1.info.short_description = "short description";
   parser1.info.url = "www.seqan.de";
   parser1.info.short_copyright = "short copyright";
   parser1.info.long_copyright = "long_copyright";
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
   EXPECT_EXIT(parser1.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");

   my_stdout = testing::internal::GetCapturedStdout();
   expected = std::string("<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01//EN\" http://www.w3.org/TR/html4/strict.dtd\">\n"
                          "<html lang=\"en\">\n"
                          "<head>\n"
                          "<meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\">\n"
                          "<title>program_full_options &mdash; short description</title>\n"
                          "</head>\n"
                          "<body>\n"
                          "<h1>program_full_options</h1>\n"
                          "<div>short description</div>\n"
                          "<h2>Synopsis</h2>\n"
                          "<p>\n"
                          "<strong>./some_binary_name</strong> synopsis\n"
                          "<br />\n"
                          "<strong>./some_binary_name</strong> synopsis2\n"
                          "<br />\n"
                          "</p>\n"
                          "<h2>Description</h2>\n"
                          "<p>\n"
                          "description\n"
                          "</p>\n"
                          "<p>\n"
                          "description2\n"
                          "</p>\n"
                          "<h2>Positional Arguments</h2>\n"
                          "<dl>\n"
                          "<dt><strong>ARGUMENT-1</strong> (<em>List</em> of <em>std::string</em>'s)</dt>\n"
                          "<dd>this is a positional option. Default: []. </dd>\n"
                          "<dt><strong>ARGUMENT-2</strong> (<em>List</em> of <em>std::string</em>'s)</dt>\n"
                          "<dd>this is a positional option. Default: []. </dd>\n"
                          "</dl>\n"
                          "<h2>Options</h2>\n"
                          "<h3>Basic options:</h3>\n"
                          "<dl>\n"
                          "<dt><strong>-h</strong>, <strong>--help</strong></dt>\n"
                          "<dd>Prints the help page.</dd>\n"
                          "<dt><strong>-hh</strong>, <strong>--advanced-help</strong></dt>\n"
                          "<dd>Prints the help page including advanced options.</dd>\n"
                          "<dt><strong>--version</strong></dt>\n"
                          "<dd>Prints the version information.</dd>\n"
                          "<dt><strong>--copyright</strong></dt>\n"
                          "<dd>Prints the copyright/license information.</dd>\n"
                          "<dt><strong>--export-help</strong> (std::string)</dt>\n"
                          "<dd>Export the help page information. Value must be one of [html, man].</dd>\n"
                          "</dl>\n"
                          "<h3></h3>\n"
                          "<dl>\n"
                          "<dt><strong>-i</strong>, <strong>--int</strong> (<em>signed 32 bit integer</em>)</dt>\n"
                          "<dd>this is a int option. Default: 5. </dd>\n"
                          "<dt><strong>-j</strong>, <strong>--jint</strong> (<em>signed 32 bit integer</em>)</dt>\n"
                          "<dd>this is a int option. Default: 5. </dd>\n"
                          "<dt><strong>-f</strong>, <strong>--flag</strong></dt>\n"
                          "<dd>this is a flag.</dd>\n"
                          "<dt><strong>-k</strong>, <strong>--kflag</strong></dt>\n"
                          "<dd>this is a flag.</dd>\n"
                          "</dl>\n"
                          "<h2>Examples</h2>\n"
                          "<p>\n"
                          "example\n"
                          "</p>\n"
                          "<p>\n"
                          "example2\n"
                          "</p>\n"
                          "<h2>Version</h2>\n"
                          "<strong>Last update:</strong> <br>\n"
                          "<strong>program_full_options version:</strong> <br>\n"
                          "<strong>SeqAn version:</strong> 3.0.0<br>\n"
                          "<h2>Url</h2>\n"
                          "www.seqan.de<br>\n"
                          "<br>\n"
                          "<h2>Legal</h2>\n"
                          "<strong>program_full_options Copyright: </strong>short copyright<br>\n"
                          "<strong>SeqAn Copyright:</strong> 2006-2019 Knut Reinert, FU-Berlin; released under the 3-clause BSDL.<br>\n"
                          "<strong>In your academic works please cite:</strong> citation<br>\n"
                          "For full copyright and/or warranty information see <tt>--copyright</tt>.\n"
                          "</body></html>");
   EXPECT_EQ(my_stdout, expected);
}

TEST(export_help, parse_error)
{
    const char * argv[]  = {"./help_add_test", "--export-help"};
    const char * argv2[] = {"./help_add_test", "--export-help=atml"};
    const char * argv3[] = {"./help_add_test", "--export-help", "atml"};

    // no value after --export-help
    EXPECT_THROW((argument_parser{"test_parser", 2, argv}), parser_invalid_argument);

    // wrong value after --export-help
    EXPECT_THROW((argument_parser{"test_parser", 2, argv2}), validation_failed);

    // wrong value after --export-help
    EXPECT_THROW((argument_parser{"test_parser", 3, argv3}), validation_failed);
}
