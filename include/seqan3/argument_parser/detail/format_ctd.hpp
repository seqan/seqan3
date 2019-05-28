// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Gianvito Urgese <gianvito.urgese AT polito.it>
 * \author Emanuele Parisi <emanuele.parisi AT polito.it>
 * \brief Contains the format_ctd struct and its helper functions.
 */

#pragma once

#include <iostream>
#include <regex>

#include <seqan3/argument_parser/detail/format_base.hpp>

namespace seqan3::detail
{

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
    
    /*!\brief Adds calls for appending 'clielement' and 'ITEM' nodes to the CTD DOM tree to be evaluated at parse time.
     *
     * \copydetails argument_parser::add_option
     */
    template<typename option_type, typename validator_type>
    void add_option(option_type & value,
                    char const short_id,
                    std::string const & long_id,
                    std::string const & desc,
                    option_spec const & spec,
                    validator_type && validator) 
    {
        // TODO (emanueleparisi) The current version of the CTD exporter does not support list options.
        if (SequenceContainer<option_type> && !std::is_same_v<option_type, std::string>)
            throw parser_design_error("At the moment, the CTD exporter does not support list options");
        
        // Do not report in the CTD file options marked as HIDDEN.
        if (spec == HIDDEN)
            return;

        // Register 'clielement' node generation callback.
        clielement_option_callbacks.push_back([this, short_id, long_id] (std::string app_name) 
        {
            if (long_id.empty())
                append_ctd_clielement_node(prepend_dash(short_id), 
                                           build_reference_name(app_name,
                                                                short_id));
            else
                append_ctd_clielement_node(prepend_dash(long_id), 
                                           build_reference_name(app_name,
                                                                long_id));
        });

        // Register ITEM callback.
        item_option_callbacks.push_back([this, value, short_id, long_id, desc, spec, validator] () 
	{
            std::string argument_name = {};
            std::string argument_required = {};
            std::string argument_advanced = {};

            // Get correct value for argument name.
            if (long_id.empty())
                argument_name = std::string{short_id};
            else
                argument_name = long_id;
	    
            // Get correct value for required/advanced flags.
            argument_required = "false";
            argument_advanced = "false";
            if (spec == REQUIRED)
            {
                argument_required = "true";
                argument_advanced = "false";
            }
            else if (spec == ADVANCED)
            {
                argument_required = "false";
                argument_advanced = "true";
            }
            
            // Register CTD 'ITEM' node creation callback.
            append_ctd_item_node(argument_name,
                                 get_type_as_ctd_string(value,
                                                        validator),
                                 desc,
                                 "",
                                 "*.*",
                                 argument_required,
                                 argument_advanced,
                                 "");
        });
    }

    /*!\brief Adds calls for appending 'clielement' and 'ITEM' nodes to the CTD DOM tree to be evaluated at parse time.
     *
     * \copydetails argument_parser::add_flag
     */
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

    /*!\brief Adds calls for appending 'clielement' and 'ITEM' nodes to the CTD DOM tree to be evaluated at parse time.
     *
     * \copydetails argument_parser::add_positional_option
     */
    template<typename option_type, typename validator_type>
    void add_positional_option(option_type & value,
                               std::string const & desc,
                               validator_type && validator) 
    {
        unsigned argument_id = args_counter;

        clielement_argument_callbacks.push_back([this, argument_id] (std::string app_name) 
	{
            append_ctd_clielement_node("",
                                       build_reference_name(app_name,
                                                            argument_id).data());
        });

        // Register ITEM callback.
        item_argument_callbacks.push_back([this, value, desc, validator, argument_id] () 
	{
            append_ctd_item_node(build_argument_name(argument_id),
                                 get_type_as_ctd_string(value,
                                                        validator),
                                 desc,
                                 "",
                                 "*.*",
                                 "true",
                                 "false",
                                 "");
        });
	
	// Increment argument counter
        args_counter++;
    }

    /*!\brief Builds the CTD document DOM tree and prints it to standard output.
     *
     * \param[in] meta The meta information that are needed for building the DOM tree.
     */
    void parse(argument_parser_meta_data const & meta) 
    {
        // Build the CTD document line by line inside the string stream.
        append_ctd_declaration_node();
        append_ctd_tool_node(meta);
        
        // Print the CTD file on the standard output stream.
        std::cout << ctd_stream.str();

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

    /*!\brief Build the identifier for the i-th positional option.
     *
     * \param[in] argument_numeric_id The numeric identifier of the argument being described.
     * \return The string representing the i-th positional option identifier.
     */
    std::string
    build_argument_name(unsigned argument_numeric_id)
    {
        return std::string{"argument-"} + std::to_string(argument_numeric_id);
    }

    /*!\brief Build the 'referenceName' attribute of the 'mapping' node in the CTD XML file.
     *
     * \param[in] app_name The name of the application the parser refers to.
     * \param[in] long_id The long identifer of the option being described.
     * \return The reference name of the option being described.
     */
    std::string 
    build_reference_name(std::string const & app_name,
                         std::string const & long_id)
    {
        return app_name + '.' + long_id;
    }

    /*!\brief Build the 'referenceName' attribute of the 'mapping' node in the CTD XML file.
     *
     * \param[in] app_name The name of the application the parser refers to.
     * \param[in] short_id The short identifer of the option being described.
     * \return The reference name of the option being described.
     */
    std::string 
    build_reference_name(std::string const & app_name,
                         char const short_id)
    {
        return app_name + '.' + short_id;
    }

    /*!\brief Build the 'referenceName' attribute of the 'mapping' node in the CTD XML file.
     *
     * \param[in] app_name The name of the application the parser refers to.
     * \param[in] argument_numeric_id The numeric identifer of the argument being described.
     * \return The reference name of the argument being described.
     */
    std::string 
    build_reference_name(std::string const & app_name,
                         unsigned argument_numeric_id)
    {
        return app_name + '.' + build_argument_name(argument_numeric_id);
    }
    
    /*!\brief Append the XML 'declaration' node to the current CTD stream buffer.
     *
     * \details 
     * The 'version' and 'encoding' attributes required by every XML declaration node are
     * hard coded to '1.0' and 'UTF-8'.
     */
    void append_ctd_declaration_node() 
    {
        ctd_stream << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    }

    /*!\brief Append a CTD 'description' node to the current CTD stream buffer.
     *
     * \param[in] meta A structure storing application meta-data.
     */
    void append_ctd_description_node(argument_parser_meta_data const & meta) 
    {
        ctd_stream << "\t<description>\n";
        if (!meta.short_description.empty())
            ctd_stream << "\t\t" << escape_special_xml_chars(meta.short_description) << "\n";
        ctd_stream << "\t</description>\n";
    }

    /*!\brief Append a CTD 'manual' node to the current CTD stream buffer.
     *
     * \param[in] meta A structure storing application meta-data.
     */
    void append_ctd_manual_node(argument_parser_meta_data const & meta) 
    {
        std::string description = {};

        // Merge description lines into a single string.
        for (auto const & l : meta.description)
            description.append(l);

        ctd_stream << "\t<manual>\n";
        if (!meta.short_description.empty())
            ctd_stream << "\t\t" << escape_special_xml_chars(description) << "\n";
        ctd_stream << "\t</manual>\n";
    }

    /*!\brief Append the CTD 'clielement' subtree to the current CTD stream buffer.
     *
     * \param[in] prefixed_option_name The option name with one or multiple dash prepended.
     * \param[in] reference_option_name The option name or argument id with the application name prepended.
     */
    void append_ctd_clielement_node(std::string const & prefixed_option_name, 
                                    std::string const & reference_option_name)
    {
        ctd_stream << "\t\t<clielement "
                   << "optionIdentifier=\"" << prefixed_option_name << "\" "
                   << "isList=\"false\""
                   << ">\n";
        ctd_stream << "\t\t\t<mapping "
                   << "referenceName=\"" << reference_option_name << "\""
                   << "/>\n";
        ctd_stream << "\t\t</clielement>\n";
    }

    /*!\brief Append the CTD 'cli' subtree to the current CTD stream buffer.
     *
     * \param[in] meta A structure storing application meta-data.
     */
    void append_ctd_cli_node(argument_parser_meta_data const & meta) 
    {
        ctd_stream << "\t<cli>\n";

	// Write command line option
        for (auto & f : clielement_option_callbacks)
            f(meta.app_name);

	// Write command line arguments
        for (auto & f : clielement_argument_callbacks)
            f(meta.app_name);
        ctd_stream << "\t</cli>\n";
    }

    /*!\brief Append the CTD 'ITEM' node to the current CTD stream buffer.
     *
     * \param[in] argument_name The 'name' attribute id of the argument being created. 
     * \param[in] argument_type The 'type' attribute of the argument being created.
     * \param[in] argument_description The 'description' attribute of the argument being created.
     * \param[in] argument_restrictions The 'restrictions' attribute of the argument being created.
     * \param[in] argument_formats The 'supported_formats' attribute of the argument being created.
     * \param[in] argument_required The 'required' attribute of the argument being created.
     * \param[in] argument_advanced The 'advanced' attribute of the argument being created.
     * \param[in] argument_value The 'value' attribute of the argument being created.
     */
    void append_ctd_item_node(std::string const & argument_name,
                              std::string const & argument_type,
                              std::string const & argument_description,
                              std::string const & argument_restrictions,
                              std::string const & argument_formats,
                              std::string const & argument_required,
                              std::string const & argument_advanced,
                              std::string const & argument_value)
    {
        ctd_stream << "\t\t\t<ITEM "
                   << "name=\"" << argument_name << "\" "
                   << "type=\"" << argument_type << "\" "
                   << "description=\"" << argument_description << "\" "
                   << "restrictions=\"" << argument_restrictions << "\" "
                   << "required=\"" << argument_required << "\" "
                   << "advanced=\"" << argument_advanced << "\" "
                   << "value=\"" << argument_value << "\"";

	// If a input/output file option is added, format requirements have to be specified
        if (argument_type == std::string{"input-file"} || argument_type == std::string{"output-file"} || 
            argument_type == std::string{"input-prefix"} || argument_type == std::string{"output-prefix"})
            ctd_stream << " supported_formats=\"" << argument_formats << "\"";

        ctd_stream << "/>\n";
    }

    /*!\brief Append the CTD 'NODE' subtree to the current CTD stream buffer.
     *
     * \param[in] meta A structure storing application meta-data.
     */
    void append_ctd_node_node(argument_parser_meta_data const & meta) 
    {
        ctd_stream << "\t\t<NODE " 
                   << "name=\"" << meta.app_name << "\" "
                   << "description=\"" << escape_special_xml_chars(meta.short_description) << "\""
                   << ">\n"; 
        
	// Run option callbacks
	for (auto f : item_option_callbacks)
            f();

	// Run argument callbacks
        for (auto f : item_argument_callbacks)
            f();

        ctd_stream << "\t\t</NODE>\n";
    }

    /*!\brief Append the CTD 'PARAMETERS' subtree to the current CTD stream buffer.
     *
     * \param[in] meta A structure storing application meta-data.
     */
    void append_ctd_parameters_node(argument_parser_meta_data const & meta) 
    {
        ctd_stream << "\t<PARAMETERS version=\"1.7.0\">\n";
        append_ctd_node_node(meta);
        ctd_stream << "\t</PARAMETERS>\n";
    }

    /*!\brief Allocate and append a CTD 'tool' subtree to the current DOM tree. 
     *
     * \param[in] meta A structure storing application meta-data.
     */
    void append_ctd_tool_node(argument_parser_meta_data const & meta) 
    {
        ctd_stream << "<tool " 
                   << "name=\"" << meta.app_name << "\" "
                   << "version=\"" << (meta.version.empty() ? "0.0.0.0" : meta.version) << "\" "
                   << "docurl=\"" << (meta.url.empty() ? "" : meta.url) << "\" "
                   << "ctdVersion=\"1.7.0\""
                   << ">\n";
        append_ctd_description_node(meta);
        append_ctd_manual_node(meta);
        append_ctd_cli_node(meta);
        append_ctd_parameters_node(meta);
        ctd_stream << "</tool>\n";
    }

    //! \brief Stream storing the CTD file lines whenever command line parsing is triggered.
    std::stringstream ctd_stream;

    //! \brief List of callbacks for generating 'clielement' subtree for options and flags.
    std::vector<std::function<void(std::string)>> clielement_option_callbacks;

    //! \brief List of callbacks for generating 'clielement' subtree for arguments.
    std::vector<std::function<void(std::string)>> clielement_argument_callbacks;

    //! \brief List of callbacks for generating 'ITEM' subtree for options and flags.
    std::vector<std::function<void(void)>> item_option_callbacks;

    //! \brief List of callbacks for generating 'ITEM' subtree for arguments.
    std::vector<std::function<void(void)>> item_argument_callbacks;

    //! \brief Number of positional arguments the user added to the argument parser.
    unsigned args_counter;
};

} // namespace seqan3
