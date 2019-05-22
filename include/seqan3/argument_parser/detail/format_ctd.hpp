// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Emanuele Parisi <emanuele.parisi AT polito.it>
 * \brief Contains the format_ctd struct and its helper functions.
 */

#pragma once

#include <iostream>
#include <regex>

#include <seqan3/argument_parser/detail/format_base.hpp>

#include <cereal/external/rapidxml/rapidxml.hpp>
#include <cereal/external/rapidxml/rapidxml_print.hpp>

namespace seqan3::detail
{

// Make a namespace alias in an anonymous namespace, such that it will not be visible from outside this file.
// TODO (emanueleparisi) Is this a good programming practice ?!?
namespace {
    //!\brief The cereal inner XML library providing facilities for handling XML documents.
    namespace rxml = cereal::rapidxml;
}

/*!\brief The format that prints the Common Tool Descriptor file to std::cout.
 * \ingroup argument_parser
 *
 * \details
 *
 * The CTD file is not written immediatelly, because the whole DOM tree
 * composing the XML document can be completely built only after the parser
 * is completely initialized. Instead, every call is stored and evaluated only
 * when format_ctd::parse() is called.
 */
class format_ctd 
{

public:

    template<typename option_type, typename validator_type>
    void add_option(option_type & value,
                    char const short_id,
                    std::string const & long_id,
                    std::string const & desc,
                    option_spec const & spec,
                    validator_type && validator) 
    {
    }

    void add_flag(bool &,
                  char const /*short_id*/,
                  std::string const & /*long_id*/,
                  std::string const & /*desc*/,
                  option_spec const & /*spec*/) 
    {
    }

    template<typename option_type, typename validator_type>
    void add_positional_option(option_type & value,
                               std::string const & desc,
                               validator_type && validator) 
    {
    }

    void parse(argument_parser_meta_data const & parser_meta) 
    {
        rxml::xml_document<> *ctd_document = new rxml::xml_document<>();

        // Validate the application name using the regex provided by the CTD XML schema.
        if (!std::regex_match(parser_meta.app_name, 
                              std::regex("[A-Za-z0-9?_\\-]+"))) {
            throw parser_design_error("CTD file requires the application name to match [A-Za-z0-9_\\-]+ regex");
        }

        // Append the XML declaration node on top of the XML document. 
        append_declaration_node(ctd_document,
                                ctd_document);

        // Append the tool subtree to the CTD document DOM tree.
        append_tool_node(ctd_document, 
                         ctd_document,
                         parser_meta);
        
        // Print the CTD file on the standard output stream.
        std::cout << *ctd_document;

        std::exit(EXIT_SUCCESS);
    }

    /*!\brief Adds an help page section to the seqan3::argument_parser.
     * \param[in] title The title of the section.
     * \details This only affects the help page and other output formats.
     */
    void add_section(std::string const &)
    {
        // Not to be implemented in CTD exporter.
    }

    /*!\brief Adds an help page subsection to the seqan3::argument_parser.
     * \param[in] title The title of the subsection.
     * \details This only affects the help page and other output formats.
     */
    void add_subsection(std::string const &)
    {
        // Not to be implemented in CTD exporter.
    }

    /*!\brief Adds an help page text line to the seqan3::argument_parser.
     * \param[in] text The text to print.
     * \param[in] line_is_paragraph Whether to insert as paragraph
     *            or just a line (Default: false).
     * \details
     * If the line is not a paragraph (false), only one line break is appended, otherwise two line breaks are appended.
     * This only affects the help page and other output formats.
     */
    void add_line(std::string const &, 
                  bool)
    {
        // Not to be implemented in CTD exporter.
    }

    /*!\brief Adds an help page list item (key-value) to the seqan3::argument_parser.
     * \param[in] key  The key of the key-value pair of the list item.
     * \param[in] desc The value of the key-value pair of the list item.
     */
    void add_list_item(std::string const &, std::string const &)
    {
        // Not to be implemented in CTD exporter.
    }

private:

    /*!\brief 
     *
     * \param[in] document
     * \param[in] parent_node
     * \param[in] version
     * \param[in] encoding
     */
    void append_declaration_node(rxml::xml_document<> *document,
                                 rxml::xml_node<> *parent_node) 
    {
        rxml::xml_node<> *declaration_node = nullptr; 

        declaration_node = document->allocate_node(rxml::node_declaration);
        declaration_node->append_attribute(document->allocate_attribute("version", 
                                                                        "1.0"));
        declaration_node->append_attribute(document->allocate_attribute("encoding", 
                                                                        "UTF-8"));
        parent_node->append_node(declaration_node);
    }

    /*!\brief 
     *
     * \param[in] document
     * \param[in] parent_node
     * \param[in] parser_meta
     */
    void append_description_node(rxml::xml_document<> *document,
                                 rxml::xml_node<> *parent_node,
                                 argument_parser_meta_data const & parser_meta) 
    {
        rxml::xml_node<> *description_node = nullptr;

        description_node = document->allocate_node(rxml::node_element, 
                                                   "description",
                                                   parser_meta.short_description.data());
        parent_node->append_node(description_node);
    }

    /*!\brief 
     *
     * \param[in] document
     * \param[in] parent_node
     * \param[in] parser_meta
     */
    void append_manual_node(rxml::xml_document<> *document, 
                            rxml::xml_node<> *parent_node, 
                            argument_parser_meta_data const & parser_meta) 
    {
        rxml::xml_node<> *manual_node = nullptr;
        std::string description = {};

        // Merge description lines into a single string.
        for (auto const & l : parser_meta.description) {
            description.append(l);
        }
        manual_node = document->allocate_node(rxml::node_element, 
                                              "manual",
                                              description.data());
        parent_node->append_node(manual_node);
    }

    /*!\brief 
     *
     * \param[in] document
     * \param[in] parent_node
     */
    void append_cli_node(rxml::xml_document<> *document,
                         rxml::xml_node<> *parent_node) 
    {
        rxml::xml_node<> *cli_node = nullptr;

        cli_node = document->allocate_node(rxml::node_element, 
                                           "cli");
        parent_node->append_node(cli_node);
    }

    /*!\brief 
     *
     * \param[in] document
     * \param[in] parent_node
     * \param[in] parser_meta
     */
    void append_node_node(rxml::xml_document<> *document, 
                          rxml::xml_node<> *parent_node, 
                          argument_parser_meta_data const & parser_meta) 
    {
        rxml::xml_node<> *node_node = nullptr;

        node_node = document->allocate_node(rxml::node_element,
                                                      "NODE");
        node_node->append_attribute(document->allocate_attribute("name", 
                                                                 parser_meta.app_name.data()));
        node_node->append_attribute(document->allocate_attribute("description",
                                                                 parser_meta.short_description.data()));
        parent_node->append_node(node_node);
    }

    /*!\brief 
     *
     * \param[in] document
     * \param[in] parent_node
     * \param[in] version
     */
    void append_parameters_node(rxml::xml_document<> *document, 
                                rxml::xml_node<> *parent_node, 
                                argument_parser_meta_data const & parser_meta) 
    {
        rxml::xml_node<> *parameters_node = nullptr;

        parameters_node = document->allocate_node(rxml::node_element, 
                                                  "PARAMETERS");
        parameters_node->append_attribute(document->allocate_attribute("version",
                                                                       "1.7.0"));
        append_node_node(document, 
                         parameters_node, 
                         parser_meta);
        parent_node->append_node(parameters_node);
    }


    /*!\brief 
     *
     * \param[in] document
     * \param[in] parent_node
     * \param[in] parser_meta
     * \param[in] ctd_version
     */
    void append_tool_node(rxml::xml_document<> *document,
                          rxml::xml_node<> *parent_node,
                          argument_parser_meta_data const & parser_meta) 
    {
        rxml::xml_node<> *tool_node = nullptr;
        
        tool_node = document->allocate_node(rxml::node_element, 
                                            "tool");
        
        // Set tool node attributes. 
        tool_node->append_attribute(document->allocate_attribute("name",
                                                                 parser_meta.app_name.data()));
        if (parser_meta.version.empty())
        {
            // App version is a mandatory attribute of the 'tool' node. If the developer does not provide any
            // data, a fake 0.0.0.0 version is used.
            tool_node->append_attribute(document->allocate_attribute("version",
                                                                     "0.0.0.0"));
        }
        else
        {
             tool_node->append_attribute(document->allocate_attribute("version",
                                                                      parser_meta.version.data()));
        }
        if (!parser_meta.url.empty())
        {
            tool_node->append_attribute(document->allocate_attribute("docurl",
                                                                     parser_meta.url.data()));
        }
        tool_node->append_attribute(document->allocate_attribute("ctdVersion",
                                                                 "1.7.0"));

        // Create and append 'description', 'manual', 'cli' and 'parameters' nodes which are children 
        // of the 'tool' node.
        append_description_node(document,
                                tool_node, 
                                parser_meta);
        append_manual_node(document,
                           tool_node, 
                           parser_meta);
        append_cli_node(document,
                        tool_node);
        append_parameters_node(document, 
                               tool_node, 
                               parser_meta);

        // Append tool node to the document DOM tree.
        parent_node->append_node(tool_node);
    }
    
    std::vector<std::function<void()>> add_option_callbacks; 
};

} // namespace seqan3
