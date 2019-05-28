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

#include <seqan3/argument_parser/validators.hpp>
#include <seqan3/argument_parser/detail/format_base.hpp>

#include <cereal/external/rapidxml/rapidxml.hpp>
#include <cereal/external/rapidxml/rapidxml_print.hpp>

namespace seqan3::detail
{

// Make a namespace alias in an anonymous namespace, such that it will not be visible from outside this file.
// TODO (emanueleparisi) Is this a good programming practice ?!?
namespace 
{
    //!\brief The cereal inner XML library providing facilities for handling XML documents.
    namespace rxml = cereal::rapidxml;
}

/*!\brief The format that prints the Common Tool Descriptor file to std::cout.
 * \ingroup argument_parser
 *
 * \details
 *
 * The CTD file is not written immediately, because the whole DOM tree
 * composing the XML document can be completely built only after the parser
 * is completely initialized. Instead, every call is stored and evaluated only
 * when format_ctd::parse() is called.
 */
class format_ctd : format_base
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
        // TODO (emanueleparisi) The current version of the CTD exporter does not support list options.
        if (SequenceContainer<option_type> && 
            !std::is_same_v<option_type,
                            std::string>)
 
        {
            throw parser_design_error("At the moment, the CTD exporter does not support list options");
        }
        
        // Do not report in the CTD file options marked as HIDDEN.
        if (spec == HIDDEN)
        {
            return;
        }

        // Register 'clielement' node generation callback.
        clielement_option_callbacks.push_back([this,
                                               short_id,
                                               long_id] (rxml::xml_document<> *pool, 
                                                         rxml::xml_node<> *parent_node,
                                                         std::string app_name) {
            char *prefixed_option_name = nullptr;
            char *reference_option_name = nullptr;
            
            // Allocate helper variables related to the DOM tree construction, getting memory 
            // from the CTD document memory pool.
            if (long_id.empty())
            {
                prefixed_option_name = pool->allocate_string(prepend_dash(short_id).data());
                reference_option_name = pool->allocate_string(prepend_app_name(app_name, 
                                                                               short_id).data());
            }
            else
            {
                prefixed_option_name = pool->allocate_string(prepend_dash(long_id).data());
                reference_option_name = pool->allocate_string(prepend_app_name(app_name,
                                                                               long_id).data());
            }

            // Build and append 'clielement' subtree.
            append_clielement_node(pool, 
                                   parent_node, 
                                   prefixed_option_name, 
                                   reference_option_name);
        });

        // Register ITEM callback.
        item_option_callbacks.push_back([this,
                                         value,
                                         short_id,
                                         long_id,
                                         desc,
                                         spec,
                                         validator] (rxml::xml_document<> *pool, 
                                                     rxml::xml_node<> *parent_node) {
            char *argument_name = nullptr;
            char *argument_type = nullptr;
            char *argument_description = nullptr;
            char *argument_restrictions = nullptr;
            char *argument_formats = nullptr;
            char *argument_required = nullptr;
            char *argument_advanced = nullptr;
            char *argument_value = nullptr;

            // Allocate helper variables related to the DOM tree construction, getting 
            // memory from the CTD document memory pool.
            if (long_id.empty())
            {
                argument_name = pool->allocate_string(std::string{short_id}.data());
            }
            else
            {
                argument_name = pool->allocate_string(long_id.data());
            }
            argument_type = pool->allocate_string(guess_gkn_type(value, 
                                                                 validator).data());
            argument_description = pool->allocate_string(desc.data()); 

            // Write the 'restrictions' and 'supported_formats' ITEM node attributes.
            // TODO (emanueleparisi) Here support for restriction and supported formats is missing ! 
            // For that to be implement we need validators to provide such information. For the sake 
            // of first releases, we ignore any constraints posed by the user.
            argument_restrictions = pool->allocate_string("");
            argument_formats = pool->allocate_string("*.*");

            // Write the 'required' and 'advanced' attributes.
            if (spec == REQUIRED)
            {
                argument_required = pool->allocate_string("true");
                argument_advanced = pool->allocate_string("false");
            }
            else if (spec == ADVANCED)
            {
                argument_required = pool->allocate_string("false");
                argument_advanced = pool->allocate_string("true");
            }
            else
            {
                argument_required = pool->allocate_string("false");
                argument_advanced = pool->allocate_string("false");
            }

            // For non-list options, append 'value' attribute to ite node. For list options, create
            // the ITEMLIST subtree. Unfortunately, the CTD exporter does not support list options.
            argument_value = pool->allocate_string("");
        
            append_item_node(pool,
                             parent_node,
                             argument_name,
                             argument_type,
                             argument_description,
                             argument_restrictions,
                             argument_formats,
                             argument_required,
                             argument_advanced,
                             argument_value);
        });
    }

    void add_flag(bool & value,
                  char const short_id,
                  std::string const & long_id,
                  std::string const & desc,
                  option_spec const & spec) 
    {
        add_option(value, 
                   short_id,
                   long_id,
                   desc,
                   spec,
                   default_validator<bool> {});
    }

    template<typename option_type, typename validator_type>
    void add_positional_option(option_type & value,
                               std::string const & desc,
                               validator_type && validator) 
    {
        clielement_argument_callbacks.push_back([this] (rxml::xml_document<> *pool, 
                                                        rxml::xml_node<> *parent_node,
                                                        argument_parser_meta_data const & meta) {
            std::string reference_option_suffix = {};
            char *reference_option_name = nullptr;

            reference_option_suffix = std::string{"argument-"}.append(std::to_string(args_counter));
            reference_option_name = pool->allocate_string(prepend_app_name(meta.app_name,
                                                                           reference_option_suffix).data());
            append_clielement_node(pool,
                                   parent_node,
                                   "",
                                   reference_option_name);
        });

        // Register ITEM callback.
        item_option_callbacks.push_back([this,
                                         value,
                                         desc,
                                         validator] (rxml::xml_document<> *pool, 
                                                     rxml::xml_node<> *parent_node) {
            char *argument_name = nullptr;
            char *argument_type = nullptr;
            char *argument_description = nullptr;
            char *argument_restrictions = nullptr;
            char *argument_formats = nullptr;
            char *argument_required = nullptr;
            char *argument_advanced = nullptr;
            char *argument_value = nullptr;
            rxml::xml_node<> *item_node = nullptr;

            // Allocate helper variables related to the DOM tree construction, getting 
            // memory from the CTD document memory pool.
            argument_name = pool->allocate_string(std::string{"argument-"}.append(std::to_string(args_counter)).data());
            argument_type = pool->allocate_string(guess_gkn_type(value, 
                                                                 validator).data());
            argument_description = pool->allocate_string(desc.data()); 

            // Write the 'restrictions' and 'supported_formats' ITEM node attributes.
            // TODO (emanueleparisi) Here support for restriction and supported formats is missing ! 
            // For that to be implement we need validators to provide such information. For the sake 
            // of first releases, we ignore any constraints posed by the user.
            item_node->append_attribute(pool->allocate_attribute("restrictions",
                                                                 ""));
            item_node->append_attribute(pool->allocate_attribute("supported_formats",
                                                                 "*.*"));

            // 'required' and 'advanced' attributes have default values for positional arguments.
            argument_required = pool->allocate_attribute("required",
                                                         "true");
            argument_required = pool->allocate_attribute("advanced",
                                                         "false");

            // For non-list options, append 'value' attribute to ite node. 
            // For list options, create the ITEMLIST subtree. Unfortunately, the CTD exporter
            // does not support list options.
            argument_value = pool->allocate_attribute("value",
                                                      "");
        
            append_item_node(pool,
                             parent_node,
                             argument_name,
                             argument_type,
                             argument_description,
                             argument_restrictions,
                             argument_formats,
                             argument_required,
                             argument_advanced,
                             argument_value);
        });

        args_counter++;
    }

    void parse(argument_parser_meta_data const & meta) 
    {
        rxml::xml_document<> *ctd_document = new rxml::xml_document<>();

        // Validate the application name using the regex provided by the CTD XML schema.
        if (!std::regex_match(meta.app_name, 
                              std::regex("[A-Za-z0-9?_\\-]+"))) {
            throw parser_design_error("CTD requires the application name to match [A-Za-z0-9_\\-]+");
        }

        // Append the XML declaration node on top of the XML document. 
        append_declaration_node(ctd_document,
                                ctd_document);

        // Append the tool subtree to the CTD document DOM tree.
        append_tool_node(ctd_document, 
                         ctd_document,
                         meta);
        
        // Print the CTD file on the standard output stream.
        std::cout << *ctd_document;

        std::exit(EXIT_SUCCESS);
    }

    // functions are not needed for command line parsing but are part of the format interface.
    //!\cond
    void add_section(std::string const &) {}
    void add_subsection(std::string const &) {}
    void add_line(std::string const &, bool) {}
    void add_list_item(std::string const &, std::string const &) {}
    //!\endcond

private:

    std::string
    prepend_app_name(std::string const & app_name,
                     std::string const & long_id)
    {
        return app_name + '.' + long_id;
    }

    std::string
    prepend_app_name(std::string const & app_name,
                     char const short_id)
    {
        return app_name + '.' + short_id;
    }

    template<typename option_type, typename validator_type>
    std::string
    guess_gkn_type(option_type const & /* option */,
                   validator_type && /* validator */)
    {
        std::string gkn_type_name = {};

        if constexpr (std::is_same_v<option_type,
                                     bool>)
        {
            gkn_type_name = "bool";
        } 
        else if constexpr (std::is_integral_v<option_type>)
        {
            gkn_type_name = "int";
        }
        else if constexpr (std::is_same_v<option_type,
                                          float>)
        {
            gkn_type_name = "float";
        }
        else if constexpr (std::is_same_v<option_type,
                                          double>)
        {
            gkn_type_name = "double";
        }
        else
        {
            if constexpr (std::is_same_v<validator_type, 
                                         input_file_validator>)
            {
                gkn_type_name = "input-file";
            }
            else if constexpr (std::is_same_v<validator_type,
                                              output_file_validator>)
            {
                gkn_type_name = "output-file";
            }
            else if constexpr (std::is_same_v<validator_type,
                                              input_directory_validator>)
            {
                gkn_type_name = "input-prefix";
            }
            else if constexpr (std::is_same_v<validator_type,
                                              output_directory_validator>)
            {
                gkn_type_name = "output-prefix";
            }
            else
            {
                gkn_type_name = "string";
            }
        }

        return gkn_type_name;
    }

    /*!\brief 
     *
     * \param[in] document
     * \param[in] parent_node
     * \param[in] version
     * \param[in] encoding
     */
    void append_declaration_node(rxml::xml_document<> *pool,
                                 rxml::xml_node<> *parent_node) 
    {
        rxml::xml_node<> *declaration_node = nullptr; 

        declaration_node = pool->allocate_node(rxml::node_declaration);
        declaration_node->append_attribute(pool->allocate_attribute("version", 
                                                                    "1.0"));
        declaration_node->append_attribute(pool->allocate_attribute("encoding", 
                                                                    "UTF-8"));
        parent_node->append_node(declaration_node);
    }

    /*!\brief 
     *
     * \param[in] document
     * \param[in] parent_node
     * \param[in] meta
     */
    void append_description_node(rxml::xml_document<> *pool,
                                 rxml::xml_node<> *parent_node,
                                 argument_parser_meta_data const & meta) 
    {
        rxml::xml_node<> *description_node = nullptr;

        description_node = pool->allocate_node(rxml::node_element, 
                                               "description",
                                               meta.short_description.data());
        parent_node->append_node(description_node);
    }

    /*!\brief 
     *
     * \param[in] document
     * \param[in] parent_node
     * \param[in] meta
     */
    void append_manual_node(rxml::xml_document<> *pool, 
                            rxml::xml_node<> *parent_node, 
                            argument_parser_meta_data const & meta) 
    {
        rxml::xml_node<> *manual_node = nullptr;
        std::string description = {};

        // Merge description lines into a single string.
        for (auto const & l : meta.description)
        {
            description.append(l);
        }
        manual_node = pool->allocate_node(rxml::node_element, 
                                          "manual",
                                          description.data());
        parent_node->append_node(manual_node);
    }

    void append_clielement_node(rxml::xml_document<> *pool, 
                                rxml::xml_node<> *parent_node, 
                                const char *prefixed_option_name, 
                                char const *reference_option_name)
    {
        rxml::xml_node<> *clielement_node = nullptr;
        rxml::xml_node<> *mapping_node = nullptr;
            
        // Allocate and fill 'clielement' node.
        clielement_node = pool->allocate_node(rxml::node_element,
                                              "clielement");
        clielement_node->append_attribute(pool->allocate_attribute("optionIdentifier",
                                                                   prefixed_option_name)); 
        
        // At the moment, list options are not supported by the CTD exporter.
        clielement_node->append_attribute(pool->allocate_attribute("isList",
                                                                   "false"));

        // Allocate and fill 'mapping' node.
        mapping_node = pool->allocate_node(rxml::node_element, 
                                           "mapping");
        mapping_node->append_attribute(pool->allocate_attribute("referenceName",
                                                                reference_option_name));

        // Build 'clielement' subtree.
        parent_node->append_node(clielement_node);
        clielement_node->append_node(mapping_node);
    }

    /*!\brief 
     *
     * \param[in] document
     * \param[in] parent_node
     * \param[in] meta
     */
    void append_cli_node(rxml::xml_document<> *pool,
                         rxml::xml_node<> *parent_node,
                         argument_parser_meta_data const & meta) 
    {
        rxml::xml_node<> *cli_node = nullptr;

        cli_node = pool->allocate_node(rxml::node_element, 
                                       "cli");
        for (auto f : clielement_option_callbacks)
        {
            f(pool,
              cli_node,
              meta.app_name);
        }
        for (auto f : clielement_argument_callbacks)
        {
            f(pool,
              cli_node,
              meta.app_name);
        }
        parent_node->append_node(cli_node);
    }

    void append_item_node(rxml::xml_document<> *pool,
                          rxml::xml_node<> *parent_node,
                          char const *argument_name,
                          char const *argument_type,
                          char const *argument_description,
                          char const *argument_restrictions,
                          char const *argument_formats,
                          char const *argument_required,
                          char const *argument_advanced,
                          char const *argument_value)
    {
        rxml::xml_node<> *item_node = nullptr;

        // Create the ITEM node for non-list options. For list options, a ITEMLIST node should
        // be created instead. Unfortunately, the CTD exporter does not support list options.
        item_node = pool->allocate_node(rxml::node_element, 
                                        "ITEM");               
        
        // Append node attributes.
        item_node->append_attribute(pool->allocate_attribute("name",
                                                             argument_name));
        item_node->append_attribute(pool->allocate_attribute("type",
                                                             argument_type));
        item_node->append_attribute(pool->allocate_attribute("description",
                                                             argument_description));
        // TODO (emanueleparisi) Here support for restriction and supported formats is missing ! 
        // For that to be implement we need validators to provide such information. For the sake 
        // of first releases, we ignore any constraints posed by the user.
        item_node->append_attribute(pool->allocate_attribute("restrictions",
                                                             argument_restrictions));
        if (argument_type == std::string{"input-file"} || 
            argument_type == std::string{"output-file"} || 
            argument_type == std::string{"input-prefix"} || 
            argument_type == std::string{"output-prefix"})
        {
            item_node->append_attribute(pool->allocate_attribute("supported_formats",
                                                                 argument_formats));
        }

        // Append 'required' and 'advanced' node attributes.
        item_node->append_attribute(pool->allocate_attribute("required",
                                                             argument_required));
        item_node->append_attribute(pool->allocate_attribute("advanced",
                                                             argument_advanced));
        item_node->append_attribute(pool->allocate_attribute("value",
                                                             argument_value));
        
        // Append the 'ITEM' subtree to the parent 'NODE' node.
        parent_node->append_node(item_node);
    }

    /*!\brief 
     *
     * \param[in] document
     * \param[in] parent_node
     * \param[in] meta
     */
    void append_node_node(rxml::xml_document<> *pool, 
                          rxml::xml_node<> *parent_node, 
                          argument_parser_meta_data const & meta) 
    {
        rxml::xml_node<> *node_node = nullptr;

        node_node = pool->allocate_node(rxml::node_element,
                                        "NODE");
        node_node->append_attribute(pool->allocate_attribute("name", 
                                                             meta.app_name.data()));
        node_node->append_attribute(pool->allocate_attribute("description",
                                                             meta.short_description.data()));
        for (auto f : item_option_callbacks)
        {
            f(pool,
              node_node);
        }
        for (auto f : item_argument_callbacks)
        {
            f(pool,
              node_node);
        }
        parent_node->append_node(node_node);
    }

    /*!\brief 
     *
     * \param[in] document
     * \param[in] parent_node
     * \param[in] version
     */
    void append_parameters_node(rxml::xml_document<> *pool, 
                                rxml::xml_node<> *parent_node, 
                                argument_parser_meta_data const & meta) 
    {
        rxml::xml_node<> *parameters_node = nullptr;

        parameters_node = pool->allocate_node(rxml::node_element, 
                                              "PARAMETERS");
        parameters_node->append_attribute(pool->allocate_attribute("version",
                                                                   "1.7.0"));
        append_node_node(pool, 
                         parameters_node, 
                         meta);
        parent_node->append_node(parameters_node);
    }

    /*!\brief 
     *
     * \param[in] document
     * \param[in] parent_node
     * \param[in] meta
     * \param[in] ctd_version
     */
    void append_tool_node(rxml::xml_document<> *pool,
                          rxml::xml_node<> *parent_node,
                          argument_parser_meta_data const & meta) 
    {
        rxml::xml_node<> *tool_node = nullptr;
        
        tool_node = pool->allocate_node(rxml::node_element, 
                                            "tool");
        
        // Set tool node attributes. 
        tool_node->append_attribute(pool->allocate_attribute("name",
                                                             meta.app_name.data()));
        if (meta.version.empty())
        {
            // App version is a mandatory attribute of the 'tool' node. If the developer 
            // does not provide any data, a fake 0.0.0.0 version is used.
            tool_node->append_attribute(pool->allocate_attribute("version",
                                                                 "0.0.0.0"));
        }
        else
        {
             tool_node->append_attribute(pool->allocate_attribute("version",
                                                                  meta.version.data()));
        }
        if (!meta.url.empty())
        {
            tool_node->append_attribute(pool->allocate_attribute("docurl",
                                                                 meta.url.data()));
        }
        tool_node->append_attribute(pool->allocate_attribute("ctdVersion",
                                                             "1.7.0"));

        // Create and append 'description', 'manual', 'cli' and 'parameters' nodes which are children 
        // of the 'tool' node.
        append_description_node(pool,
                                tool_node, 
                                meta);
        append_manual_node(pool,
                           tool_node, 
                           meta);
        append_cli_node(pool,
                        tool_node,
                        meta);
        append_parameters_node(pool, 
                               tool_node, 
                               meta);

        // Append tool node to the document DOM tree.
        parent_node->append_node(tool_node);
    }
   
    std::vector<std::function<void(rxml::xml_document<>*, 
                                   rxml::xml_node<>*,
                                   std::string)>> clielement_option_callbacks;
    std::vector<std::function<void(rxml::xml_document<>*, 
                                   rxml::xml_node<>*,
                                   std::string)>> clielement_argument_callbacks;
    std::vector<std::function<void(rxml::xml_document<>*, 
                                   rxml::xml_node<>*)>> item_option_callbacks;
    std::vector<std::function<void(rxml::xml_document<>*, 
                                   rxml::xml_node<>*)>> item_argument_callbacks;
    unsigned args_counter;
};

} // namespace seqan3
